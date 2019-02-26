% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:47:36
% EndTime: 2019-02-26 20:47:36
% DurationCPUTime: 0.12s
% Computational Cost: add. (354->25), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->41)
t55 = qJ(3) + pkin(9);
t54 = cos(t55);
t75 = t54 ^ 2;
t53 = sin(t55);
t59 = sin(qJ(5));
t62 = cos(qJ(1));
t65 = t62 * t59;
t60 = sin(qJ(1));
t61 = cos(qJ(5));
t66 = t60 * t61;
t47 = t53 * t65 + t66;
t69 = t54 * t59;
t42 = atan2(t47, t69);
t38 = sin(t42);
t39 = cos(t42);
t37 = t38 * t47 + t39 * t69;
t36 = 0.1e1 / t37 ^ 2;
t64 = t62 * t61;
t67 = t60 * t59;
t45 = t53 * t67 - t64;
t74 = t36 * t45;
t72 = t39 * t47;
t71 = t45 ^ 2 * t36;
t51 = 0.1e1 / t54;
t56 = 0.1e1 / t59;
t70 = t51 * t56;
t68 = t54 * t60;
t46 = t53 * t66 + t65;
t44 = 0.1e1 / t46 ^ 2;
t63 = t60 ^ 2 * t75 * t44;
t57 = 0.1e1 / t59 ^ 2;
t52 = 0.1e1 / t75;
t48 = t53 * t64 - t67;
t43 = 0.1e1 / t46;
t41 = 0.1e1 / (t47 ^ 2 * t52 * t57 + 0.1e1);
t40 = 0.1e1 / (0.1e1 + t63);
t35 = 0.1e1 / t37;
t34 = (t47 * t52 * t53 * t56 + t62) * t41;
t33 = 0.1e1 / (0.1e1 + t71);
t32 = (-t47 * t57 * t61 + t48 * t56) * t51 * t41;
t1 = [-t45 * t41 * t70, 0, t34, 0, t32, 0; (t47 * t35 - (-t38 + (-t70 * t72 + t38) * t41) * t71) * t33, 0 (-t34 * t72 * t74 + (t35 * t68 - (-t39 * t53 + (-t34 + t62) * t54 * t38) * t74) * t59) * t33, 0 (t46 * t35 - (t39 * t54 * t61 + t38 * t48 + (-t38 * t69 + t72) * t32) * t74) * t33, 0; (-t44 * t48 * t60 + t43 * t62) * t54 * t40, 0 (-t43 * t53 * t60 - t61 * t63) * t40, 0, t45 * t44 * t40 * t68, 0;];
Ja_rot  = t1;
