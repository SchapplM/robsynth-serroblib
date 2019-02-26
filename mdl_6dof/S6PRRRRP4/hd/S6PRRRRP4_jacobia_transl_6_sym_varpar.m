% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRP4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:53
% EndTime: 2019-02-26 20:16:54
% DurationCPUTime: 0.22s
% Computational Cost: add. (343->53), mult. (630->93), div. (0->0), fcn. (805->12), ass. (0->43)
t69 = pkin(5) + r_i_i_C(1);
t61 = r_i_i_C(3) + qJ(6);
t45 = cos(qJ(4));
t36 = t45 * pkin(4) + pkin(3);
t43 = sin(qJ(3));
t46 = cos(qJ(3));
t68 = r_i_i_C(2) + pkin(10) + pkin(9);
t70 = t36 * t46 + t68 * t43 + pkin(2);
t39 = qJ(4) + qJ(5);
t37 = sin(t39);
t67 = t37 * t46;
t38 = cos(t39);
t66 = t38 * t46;
t41 = sin(pkin(6));
t65 = t41 * t43;
t64 = t41 * t46;
t47 = cos(qJ(2));
t63 = t41 * t47;
t62 = t46 * t47;
t60 = cos(pkin(6));
t59 = cos(pkin(11));
t42 = sin(qJ(4));
t58 = pkin(4) * t42 + pkin(8);
t40 = sin(pkin(11));
t56 = t40 * t60;
t55 = t41 * t59;
t44 = sin(qJ(2));
t52 = t60 * t59;
t29 = t40 * t47 + t44 * t52;
t21 = t29 * t46 - t43 * t55;
t28 = t40 * t44 - t47 * t52;
t7 = t21 * t37 - t28 * t38;
t54 = t61 * (t21 * t38 + t28 * t37) - t69 * t7;
t31 = -t44 * t56 + t59 * t47;
t23 = t31 * t46 + t40 * t65;
t30 = t59 * t44 + t47 * t56;
t9 = t23 * t37 - t30 * t38;
t53 = -t69 * t9 + t61 * (t23 * t38 + t30 * t37);
t33 = t60 * t43 + t44 * t64;
t18 = t33 * t37 + t38 * t63;
t51 = t61 * (t33 * t38 - t37 * t63) - t69 * t18;
t49 = t61 * t37 + t69 * t38 + t36;
t1 = [0, t58 * t31 + t69 * (-t30 * t66 + t31 * t37) + t61 * (-t30 * t67 - t31 * t38) - t70 * t30, t68 * t23 + t49 * (-t31 * t43 + t40 * t64) (-t23 * t42 + t30 * t45) * pkin(4) + t53, t53, t9; 0, t58 * t29 + t69 * (-t28 * t66 + t29 * t37) + t61 * (-t28 * t67 - t29 * t38) - t70 * t28, t68 * t21 + t49 * (-t29 * t43 - t46 * t55) (-t21 * t42 + t28 * t45) * pkin(4) + t54, t54, t7; 1 (t69 * (t37 * t44 + t38 * t62) + t61 * (t37 * t62 - t38 * t44) + t58 * t44 + t70 * t47) * t41, t68 * t33 + t49 * (-t44 * t65 + t60 * t46) (-t33 * t42 - t45 * t63) * pkin(4) + t51, t51, t18;];
Ja_transl  = t1;
