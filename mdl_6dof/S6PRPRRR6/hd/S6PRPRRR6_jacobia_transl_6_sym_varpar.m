% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:55
% EndTime: 2019-02-26 19:56:55
% DurationCPUTime: 0.13s
% Computational Cost: add. (205->43), mult. (401->74), div. (0->0), fcn. (506->12), ass. (0->34)
t23 = sin(pkin(11));
t25 = cos(pkin(11));
t32 = cos(qJ(2));
t26 = cos(pkin(6));
t29 = sin(qJ(2));
t38 = t26 * t29;
t15 = -t23 * t38 + t25 * t32;
t22 = qJ(5) + qJ(6);
t20 = sin(t22);
t21 = cos(t22);
t37 = t26 * t32;
t14 = t23 * t37 + t25 * t29;
t28 = sin(qJ(4));
t24 = sin(pkin(6));
t31 = cos(qJ(4));
t40 = t24 * t31;
t8 = t14 * t28 + t23 * t40;
t46 = (t15 * t21 - t8 * t20) * r_i_i_C(1) + (-t15 * t20 - t8 * t21) * r_i_i_C(2);
t12 = t23 * t29 - t25 * t37;
t10 = -t12 * t28 + t25 * t40;
t13 = t23 * t32 + t25 * t38;
t45 = (t10 * t20 + t13 * t21) * r_i_i_C(1) + (t10 * t21 - t13 * t20) * r_i_i_C(2);
t39 = t24 * t32;
t17 = t26 * t31 - t28 * t39;
t41 = t24 * t29;
t44 = (-t17 * t20 + t21 * t41) * r_i_i_C(1) + (-t17 * t21 - t20 * t41) * r_i_i_C(2);
t43 = r_i_i_C(3) + pkin(10) + pkin(9);
t42 = t24 * t28;
t30 = cos(qJ(5));
t36 = t30 * pkin(5) + r_i_i_C(1) * t21 - r_i_i_C(2) * t20 + pkin(4);
t27 = sin(qJ(5));
t35 = t27 * pkin(5) + t20 * r_i_i_C(1) + t21 * r_i_i_C(2) + pkin(2) + pkin(8);
t34 = t36 * t28 - t43 * t31 + qJ(3);
t1 = [0, -t35 * t14 + t34 * t15, t14, t43 * t8 + t36 * (t14 * t31 - t23 * t42) (t15 * t30 - t27 * t8) * pkin(5) + t46, t46; 0, -t35 * t12 + t34 * t13, t12, -t43 * t10 + t36 * (t12 * t31 + t25 * t42) (t10 * t27 + t13 * t30) * pkin(5) + t45, t45; 1 (t34 * t29 + t35 * t32) * t24, -t39, t43 * t17 + t36 * (-t26 * t28 - t31 * t39) (-t17 * t27 + t30 * t41) * pkin(5) + t44, t44;];
Ja_transl  = t1;
