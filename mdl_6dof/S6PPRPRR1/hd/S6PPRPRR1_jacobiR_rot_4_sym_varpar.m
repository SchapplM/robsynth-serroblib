% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRPRR1_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobiR_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:48
% EndTime: 2019-02-26 19:39:48
% DurationCPUTime: 0.07s
% Computational Cost: add. (38->18), mult. (110->40), div. (0->0), fcn. (154->12), ass. (0->27)
t75 = sin(pkin(11));
t77 = sin(pkin(6));
t88 = t75 * t77;
t82 = cos(pkin(6));
t87 = t75 * t82;
t80 = cos(pkin(11));
t86 = t77 * t80;
t85 = t80 * t82;
t73 = sin(pkin(13));
t78 = cos(pkin(13));
t83 = sin(qJ(3));
t84 = cos(qJ(3));
t72 = -t84 * t73 - t83 * t78;
t71 = t83 * t73 - t84 * t78;
t81 = cos(pkin(7));
t79 = cos(pkin(12));
t76 = sin(pkin(7));
t74 = sin(pkin(12));
t70 = -t74 * t87 + t80 * t79;
t69 = -t80 * t74 - t79 * t87;
t68 = t74 * t85 + t75 * t79;
t67 = -t75 * t74 + t79 * t85;
t66 = t72 * t81;
t65 = t71 * t81;
t64 = t72 * t76;
t63 = t71 * t76;
t1 = [0, 0, -t63 * t88 - t69 * t65 + t70 * t72, 0, 0, 0; 0, 0, t63 * t86 - t67 * t65 + t68 * t72, 0, 0, 0; 0, 0, -t82 * t63 + (-t65 * t79 + t72 * t74) * t77, 0, 0, 0; 0, 0, t64 * t88 + t69 * t66 + t70 * t71, 0, 0, 0; 0, 0, -t64 * t86 + t67 * t66 + t68 * t71, 0, 0, 0; 0, 0, t82 * t64 + (t66 * t79 + t71 * t74) * t77, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
