% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:43
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPP8_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_jacobiR_rot_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:42:59
% EndTime: 2019-02-22 10:42:59
% DurationCPUTime: 0.06s
% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
t75 = sin(qJ(3));
t76 = sin(qJ(1));
t86 = t76 * t75;
t77 = cos(qJ(4));
t85 = t76 * t77;
t74 = sin(qJ(4));
t78 = cos(qJ(3));
t84 = t78 * t74;
t83 = t78 * t77;
t79 = cos(qJ(1));
t82 = t79 * t75;
t81 = t79 * t77;
t80 = t79 * t78;
t73 = t76 * t74 - t75 * t81;
t72 = t74 * t82 + t85;
t71 = t79 * t74 + t75 * t85;
t70 = t74 * t86 - t81;
t1 = [-t80, 0, t86, 0, 0, 0; -t76 * t78, 0, -t82, 0, 0, 0; 0, 0, t78, 0, 0, 0; t73, 0, -t76 * t83, t70, 0, 0; -t71, 0, t77 * t80, -t72, 0, 0; 0, 0, t75 * t77, t84, 0, 0; t72, 0, t76 * t84, t71, 0, 0; t70, 0, -t74 * t80, t73, 0, 0; 0, 0, -t75 * t74, t83, 0, 0;];
JR_rot  = t1;
