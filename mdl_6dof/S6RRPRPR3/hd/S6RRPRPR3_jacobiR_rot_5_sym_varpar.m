% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:24
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:23:52
% EndTime: 2019-02-22 11:23:52
% DurationCPUTime: 0.04s
% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t82 = qJ(2) + pkin(10);
t78 = sin(t82);
t83 = sin(qJ(1));
t90 = t83 * t78;
t81 = qJ(4) + pkin(11);
t79 = cos(t81);
t89 = t83 * t79;
t80 = cos(t82);
t88 = t83 * t80;
t84 = cos(qJ(1));
t87 = t84 * t78;
t86 = t84 * t79;
t85 = t84 * t80;
t77 = sin(t81);
t76 = t83 * t77 + t79 * t85;
t75 = -t77 * t85 + t89;
t74 = t84 * t77 - t79 * t88;
t73 = t77 * t88 + t86;
t1 = [t74, -t78 * t86, 0, t75, 0, 0; t76, -t78 * t89, 0, -t73, 0, 0; 0, t80 * t79, 0, -t78 * t77, 0, 0; t73, t77 * t87, 0, -t76, 0, 0; t75, t77 * t90, 0, t74, 0, 0; 0, -t80 * t77, 0, -t78 * t79, 0, 0; -t90, t85, 0, 0, 0, 0; t87, t88, 0, 0, 0, 0; 0, t78, 0, 0, 0, 0;];
JR_rot  = t1;
