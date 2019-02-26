% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR2
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
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:11
% EndTime: 2019-02-26 19:54:12
% DurationCPUTime: 0.21s
% Computational Cost: add. (328->49), mult. (723->86), div. (0->0), fcn. (949->14), ass. (0->39)
t32 = sin(pkin(12));
t39 = sin(qJ(2));
t42 = cos(qJ(2));
t51 = cos(pkin(12));
t47 = -t39 * t32 + t42 * t51;
t38 = sin(qJ(4));
t41 = cos(qJ(4));
t31 = qJ(5) + qJ(6);
t29 = sin(t31);
t30 = cos(t31);
t40 = cos(qJ(5));
t48 = t40 * pkin(5) + r_i_i_C(1) * t30 - r_i_i_C(2) * t29 + pkin(4);
t56 = r_i_i_C(3) + pkin(10) + pkin(9);
t44 = t56 * t38 + t48 * t41 + pkin(3);
t25 = -t42 * t32 - t39 * t51;
t33 = sin(pkin(11));
t35 = cos(pkin(11));
t36 = cos(pkin(6));
t45 = t47 * t36;
t12 = t33 * t25 + t35 * t45;
t23 = t25 * t36;
t13 = -t35 * t23 + t33 * t47;
t34 = sin(pkin(6));
t54 = t35 * t34;
t8 = t13 * t41 - t38 * t54;
t59 = (-t12 * t30 - t8 * t29) * r_i_i_C(1) + (t12 * t29 - t8 * t30) * r_i_i_C(2);
t49 = -t33 * t23 - t35 * t47;
t55 = t33 * t34;
t10 = t38 * t55 - t41 * t49;
t15 = t35 * t25 - t33 * t45;
t58 = (-t10 * t29 - t15 * t30) * r_i_i_C(1) + (-t10 * t30 + t15 * t29) * r_i_i_C(2);
t22 = t25 * t34;
t18 = -t22 * t41 + t36 * t38;
t21 = t47 * t34;
t57 = (-t18 * t29 - t21 * t30) * r_i_i_C(1) + (-t18 * t30 + t21 * t29) * r_i_i_C(2);
t53 = t36 * t42;
t37 = sin(qJ(5));
t46 = t37 * pkin(5) + t29 * r_i_i_C(1) + t30 * r_i_i_C(2) + pkin(8);
t1 = [0 (-t33 * t53 - t35 * t39) * pkin(2) - t46 * t49 + t44 * t15, t55, t56 * t10 + t48 * (t38 * t49 + t41 * t55) (-t10 * t37 - t15 * t40) * pkin(5) + t58, t58; 0 (-t33 * t39 + t35 * t53) * pkin(2) + t46 * t13 + t44 * t12, -t54, t56 * t8 + t48 * (-t13 * t38 - t41 * t54) (-t12 * t40 - t37 * t8) * pkin(5) + t59, t59; 1, t34 * t42 * pkin(2) + t44 * t21 - t46 * t22, t36, t56 * t18 + t48 * (t22 * t38 + t36 * t41) (-t18 * t37 - t21 * t40) * pkin(5) + t57, t57;];
Ja_transl  = t1;
