% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:04:41
% EndTime: 2019-02-26 20:04:42
% DurationCPUTime: 0.20s
% Computational Cost: add. (289->51), mult. (415->85), div. (0->0), fcn. (520->14), ass. (0->38)
t24 = qJ(3) + pkin(12);
t20 = sin(t24);
t21 = cos(t24);
t35 = cos(qJ(3));
t25 = qJ(5) + qJ(6);
t22 = sin(t25);
t23 = cos(t25);
t34 = cos(qJ(5));
t40 = t34 * pkin(5) + r_i_i_C(1) * t23 - r_i_i_C(2) * t22 + pkin(4);
t48 = r_i_i_C(3) + pkin(10) + pkin(9);
t52 = t35 * pkin(3) + t48 * t20 + t40 * t21 + pkin(2);
t26 = sin(pkin(11));
t28 = cos(pkin(11));
t33 = sin(qJ(2));
t29 = cos(pkin(6));
t36 = cos(qJ(2));
t41 = t29 * t36;
t13 = t26 * t33 - t28 * t41;
t42 = t29 * t33;
t14 = t26 * t36 + t28 * t42;
t27 = sin(pkin(6));
t46 = t27 * t28;
t8 = t14 * t21 - t20 * t46;
t51 = (t13 * t23 - t8 * t22) * r_i_i_C(1) + (-t13 * t22 - t8 * t23) * r_i_i_C(2);
t16 = -t26 * t42 + t28 * t36;
t47 = t26 * t27;
t10 = t16 * t21 + t20 * t47;
t15 = t26 * t41 + t28 * t33;
t50 = (-t10 * t22 + t15 * t23) * r_i_i_C(1) + (-t10 * t23 - t15 * t22) * r_i_i_C(2);
t45 = t27 * t33;
t12 = t29 * t20 + t21 * t45;
t43 = t27 * t36;
t49 = (-t12 * t22 - t23 * t43) * r_i_i_C(1) + (-t12 * t23 + t22 * t43) * r_i_i_C(2);
t44 = t27 * t35;
t31 = sin(qJ(5));
t39 = t31 * pkin(5) + t22 * r_i_i_C(1) + t23 * r_i_i_C(2) + pkin(8) + qJ(4);
t32 = sin(qJ(3));
t1 = [0, -t15 * t52 + t39 * t16, t48 * t10 + (-t16 * t32 + t26 * t44) * pkin(3) + t40 * (-t16 * t20 + t21 * t47) t15 (-t10 * t31 + t15 * t34) * pkin(5) + t50, t50; 0, -t13 * t52 + t39 * t14, t48 * t8 + (-t14 * t32 - t28 * t44) * pkin(3) + t40 * (-t14 * t20 - t21 * t46) t13 (t13 * t34 - t31 * t8) * pkin(5) + t51, t51; 1 (t39 * t33 + t52 * t36) * t27, t48 * t12 + (t29 * t35 - t32 * t45) * pkin(3) + t40 * (-t20 * t45 + t29 * t21) -t43 (-t12 * t31 - t34 * t43) * pkin(5) + t49, t49;];
Ja_transl  = t1;
