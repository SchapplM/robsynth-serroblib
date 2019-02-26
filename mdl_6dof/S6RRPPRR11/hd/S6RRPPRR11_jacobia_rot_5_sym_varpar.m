% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR11
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR11_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:19
% EndTime: 2019-02-26 21:34:19
% DurationCPUTime: 0.09s
% Computational Cost: add. (197->25), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->41)
t59 = sin(qJ(2));
t61 = cos(qJ(2));
t62 = cos(qJ(1));
t60 = sin(qJ(1));
t66 = cos(pkin(6));
t64 = t60 * t66;
t48 = t62 * t59 + t61 * t64;
t55 = pkin(11) + qJ(5);
t52 = sin(t55);
t53 = cos(t55);
t58 = sin(pkin(6));
t68 = t58 * t60;
t40 = t48 * t52 + t53 * t68;
t38 = 0.1e1 / t40 ^ 2;
t39 = -t48 * t53 + t52 * t68;
t73 = t38 * t39;
t63 = t62 * t66;
t46 = t59 * t63 + t60 * t61;
t69 = t58 * t59;
t44 = atan2(-t46, t69);
t42 = cos(t44);
t72 = t42 * t46;
t41 = sin(t44);
t35 = -t41 * t46 + t42 * t69;
t34 = 0.1e1 / t35 ^ 2;
t49 = -t59 * t64 + t62 * t61;
t71 = t49 ^ 2 * t34;
t54 = 0.1e1 / t58;
t56 = 0.1e1 / t59;
t70 = t54 * t56;
t67 = t58 * t62;
t65 = t39 ^ 2 * t38 + 0.1e1;
t57 = 0.1e1 / t59 ^ 2;
t45 = t60 * t59 - t61 * t63;
t43 = 0.1e1 / (0.1e1 + t46 ^ 2 / t58 ^ 2 * t57);
t37 = 0.1e1 / t40;
t36 = 0.1e1 / t65;
t33 = 0.1e1 / t35;
t32 = 0.1e1 / (0.1e1 + t71);
t31 = (t46 * t57 * t61 + t45 * t56) * t54 * t43;
t1 = [-t49 * t43 * t70, t31, 0, 0, 0, 0; (-t46 * t33 - (-t41 + (t70 * t72 + t41) * t43) * t71) * t32 (-t48 * t33 - (t42 * t58 * t61 + t41 * t45 + (-t41 * t69 - t72) * t31) * t49 * t34) * t32, 0, 0, 0, 0; ((t45 * t53 + t52 * t67) * t37 - (-t45 * t52 + t53 * t67) * t73) * t36 (-t53 * t37 - t52 * t73) * t49 * t36, 0, 0, t65 * t36, 0;];
Ja_rot  = t1;
