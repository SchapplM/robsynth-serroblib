% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR7_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:57:34
% EndTime: 2019-02-26 21:57:35
% DurationCPUTime: 0.17s
% Computational Cost: add. (307->23), mult. (876->55), div. (85->9), fcn. (1282->11), ass. (0->37)
t68 = sin(qJ(4));
t71 = cos(qJ(2));
t87 = sin(qJ(2));
t88 = cos(qJ(4));
t61 = -t71 * t68 + t87 * t88;
t69 = sin(qJ(1));
t53 = t61 * t69;
t60 = t87 * t68 + t71 * t88;
t48 = atan2(t53, t60);
t45 = sin(t48);
t46 = cos(t48);
t43 = t45 * t53 + t46 * t60;
t42 = 0.1e1 / t43 ^ 2;
t89 = cos(qJ(1));
t57 = t61 * t89;
t81 = t57 ^ 2 * t42;
t40 = 0.1e1 / (0.1e1 + t81);
t41 = 0.1e1 / t43;
t55 = t60 * t69;
t56 = t60 * t89;
t84 = t46 * t53;
t59 = 0.1e1 / t60 ^ 2;
t47 = 0.1e1 / (t53 ^ 2 * t59 + 0.1e1);
t58 = 0.1e1 / t60;
t91 = (-t53 * t59 * t61 - t55 * t58) * t47;
t94 = ((-t45 * t55 + t46 * t61 + (-t45 * t60 + t84) * t91) * t42 * t57 + t56 * t41) * t40;
t67 = sin(qJ(5));
t70 = cos(qJ(5));
t52 = t56 * t70 - t69 * t67;
t50 = 0.1e1 / t52 ^ 2;
t51 = t56 * t67 + t69 * t70;
t79 = t51 ^ 2 * t50 + 0.1e1;
t44 = 0.1e1 / t79;
t49 = 0.1e1 / t52;
t83 = t50 * t51;
t90 = (-t67 * t49 + t70 * t83) * t44 * t57;
t1 = [t57 * t58 * t47, -t91, 0, t91, 0, 0; (t53 * t41 - (-t45 + (-t58 * t84 + t45) * t47) * t81) * t40, -t94, 0, t94, 0, 0; ((-t55 * t67 + t89 * t70) * t49 - (-t55 * t70 - t89 * t67) * t83) * t44, t90, 0, -t90, t79 * t44, 0;];
Ja_rot  = t1;
