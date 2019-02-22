% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:08
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPPR5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:08:34
% EndTime: 2019-02-22 11:08:34
% DurationCPUTime: 0.05s
% Computational Cost: add. (33->18), mult. (108->32), div. (0->0), fcn. (161->8), ass. (0->23)
t80 = sin(qJ(2));
t77 = sin(pkin(9));
t78 = cos(pkin(9));
t79 = sin(qJ(6));
t82 = cos(qJ(6));
t87 = t77 * t79 - t78 * t82;
t86 = t80 * t87;
t81 = sin(qJ(1));
t83 = cos(qJ(2));
t92 = t81 * t83;
t84 = cos(qJ(1));
t91 = t84 * t83;
t72 = t77 * t92 + t84 * t78;
t73 = -t84 * t77 + t78 * t92;
t90 = -t72 * t82 - t73 * t79;
t89 = t72 * t79 - t73 * t82;
t88 = t77 * t82 + t78 * t79;
t85 = t88 * t80;
t75 = t81 * t77 + t78 * t91;
t74 = t77 * t91 - t81 * t78;
t71 = t74 * t82 + t75 * t79;
t70 = -t74 * t79 + t75 * t82;
t1 = [t90, -t84 * t85, 0, 0, 0, t70; t71, -t81 * t85, 0, 0, 0, -t89; 0, t88 * t83, 0, 0, 0, -t86; t89, t84 * t86, 0, 0, 0, -t71; t70, t81 * t86, 0, 0, 0, t90; 0, -t87 * t83, 0, 0, 0, -t85; -t81 * t80, t91, 0, 0, 0, 0; t84 * t80, t92, 0, 0, 0, 0; 0, t80, 0, 0, 0, 0;];
JR_rot  = t1;
