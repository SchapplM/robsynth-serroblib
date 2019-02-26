% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:43
% EndTime: 2019-02-26 19:53:43
% DurationCPUTime: 0.20s
% Computational Cost: add. (366->49), mult. (715->85), div. (0->0), fcn. (934->14), ass. (0->39)
t67 = pkin(10) + r_i_i_C(3);
t41 = sin(qJ(6));
t44 = cos(qJ(6));
t66 = t44 * r_i_i_C(1) - t41 * r_i_i_C(2) + pkin(5);
t36 = sin(pkin(12));
t43 = sin(qJ(2));
t46 = cos(qJ(2));
t58 = cos(pkin(12));
t53 = -t43 * t36 + t46 * t58;
t35 = qJ(4) + qJ(5);
t33 = sin(t35);
t34 = cos(t35);
t45 = cos(qJ(4));
t48 = t45 * pkin(4) + t67 * t33 + t66 * t34 + pkin(3);
t37 = sin(pkin(11));
t38 = sin(pkin(6));
t62 = t37 * t38;
t39 = cos(pkin(11));
t61 = t39 * t38;
t40 = cos(pkin(6));
t60 = t40 * t46;
t28 = -t46 * t36 - t43 * t58;
t26 = t28 * t40;
t13 = -t39 * t26 + t37 * t53;
t56 = -t37 * t26 - t39 * t53;
t54 = t41 * r_i_i_C(1) + t44 * r_i_i_C(2) + pkin(8) + pkin(9);
t8 = t13 * t34 - t33 * t61;
t52 = t67 * t8 + t66 * (-t13 * t33 - t34 * t61);
t10 = t33 * t62 - t34 * t56;
t51 = t67 * t10 + t66 * (t33 * t56 + t34 * t62);
t50 = t53 * t40;
t25 = t28 * t38;
t21 = -t25 * t34 + t40 * t33;
t49 = t67 * t21 + t66 * (t25 * t33 + t40 * t34);
t42 = sin(qJ(4));
t24 = t53 * t38;
t15 = t39 * t28 - t37 * t50;
t12 = t37 * t28 + t39 * t50;
t1 = [0 (-t37 * t60 - t39 * t43) * pkin(2) - t54 * t56 + t48 * t15, t62 (t42 * t56 + t45 * t62) * pkin(4) + t51, t51 (-t10 * t41 - t15 * t44) * r_i_i_C(1) + (-t10 * t44 + t15 * t41) * r_i_i_C(2); 0 (-t37 * t43 + t39 * t60) * pkin(2) + t54 * t13 + t48 * t12, -t61 (-t13 * t42 - t45 * t61) * pkin(4) + t52, t52 (-t12 * t44 - t8 * t41) * r_i_i_C(1) + (t12 * t41 - t8 * t44) * r_i_i_C(2); 1, t38 * t46 * pkin(2) + t48 * t24 - t54 * t25, t40 (t25 * t42 + t40 * t45) * pkin(4) + t49, t49 (-t21 * t41 - t24 * t44) * r_i_i_C(1) + (-t21 * t44 + t24 * t41) * r_i_i_C(2);];
Ja_transl  = t1;
