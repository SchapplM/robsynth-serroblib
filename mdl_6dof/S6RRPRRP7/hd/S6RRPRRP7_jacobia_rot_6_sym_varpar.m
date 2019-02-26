% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:49:25
% EndTime: 2019-02-26 21:49:25
% DurationCPUTime: 0.27s
% Computational Cost: add. (491->30), mult. (1398->74), div. (142->11), fcn. (2036->11), ass. (0->45)
t109 = cos(qJ(1));
t107 = sin(qJ(2));
t108 = cos(qJ(4));
t86 = sin(qJ(4));
t89 = cos(qJ(2));
t76 = t107 * t86 + t89 * t108;
t87 = sin(qJ(1));
t70 = t76 * t87;
t85 = sin(qJ(5));
t88 = cos(qJ(5));
t63 = -t109 * t88 + t70 * t85;
t77 = t107 * t108 - t89 * t86;
t97 = t77 * t85;
t60 = atan2(-t63, t97);
t58 = cos(t60);
t104 = t58 * t63;
t57 = sin(t60);
t55 = -t57 * t63 + t58 * t97;
t54 = 0.1e1 / t55 ^ 2;
t71 = t76 * t109;
t66 = t71 * t85 + t87 * t88;
t106 = t54 * t66;
t75 = 0.1e1 / t77 ^ 2;
t84 = 0.1e1 / t85 ^ 2;
t59 = 0.1e1 / (t63 ^ 2 * t75 * t84 + 0.1e1);
t69 = t77 * t87;
t74 = 0.1e1 / t77;
t83 = 0.1e1 / t85;
t110 = (t63 * t75 * t76 * t83 + t69 * t74) * t59;
t101 = t66 ^ 2 * t54;
t52 = 0.1e1 / (0.1e1 + t101);
t53 = 0.1e1 / t55;
t72 = t77 * t109;
t117 = (t110 * t104 * t106 + ((-t58 * t76 + (t110 * t77 - t69) * t57) * t106 - t72 * t53) * t85) * t52;
t67 = t71 * t88 - t87 * t85;
t62 = 0.1e1 / t67 ^ 2;
t100 = t72 ^ 2 * t62;
t56 = 0.1e1 / (0.1e1 + t100);
t61 = 0.1e1 / t67;
t115 = (t88 * t100 + t61 * t71) * t56;
t102 = t62 * t72;
t98 = t74 * t83;
t65 = t109 * t85 + t70 * t88;
t49 = (t63 * t84 * t88 - t65 * t83) * t74 * t59;
t1 = [-t66 * t59 * t98, t110, 0, -t110, t49, 0; (-t63 * t53 - (-t57 + (t98 * t104 + t57) * t59) * t101) * t52, t117, 0, -t117 (t67 * t53 - (t58 * t77 * t88 - t57 * t65 + (-t57 * t97 - t104) * t49) * t106) * t52, 0; (t65 * t102 - t69 * t61) * t56, t115, 0, -t115, t66 * t56 * t102, 0;];
Ja_rot  = t1;
