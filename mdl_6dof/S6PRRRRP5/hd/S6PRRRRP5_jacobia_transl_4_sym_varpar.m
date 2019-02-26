% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRP5_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP5_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_jacobia_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:17:28
% EndTime: 2019-02-26 20:17:28
% DurationCPUTime: 0.16s
% Computational Cost: add. (168->57), mult. (471->110), div. (0->0), fcn. (606->12), ass. (0->43)
t50 = r_i_i_C(3) + pkin(10);
t22 = sin(pkin(7));
t49 = t22 * pkin(9);
t23 = sin(pkin(6));
t48 = t22 * t23;
t26 = cos(pkin(6));
t47 = t22 * t26;
t27 = sin(qJ(4));
t46 = t22 * t27;
t30 = cos(qJ(4));
t45 = t22 * t30;
t25 = cos(pkin(7));
t44 = t23 * t25;
t28 = sin(qJ(3));
t43 = t25 * t28;
t31 = cos(qJ(3));
t42 = t25 * t31;
t29 = sin(qJ(2));
t41 = t26 * t29;
t32 = cos(qJ(2));
t40 = t26 * t32;
t39 = t28 * t29;
t38 = t28 * t32;
t37 = t29 * t31;
t36 = t31 * t32;
t35 = t30 * r_i_i_C(1) - t27 * r_i_i_C(2) + pkin(3);
t21 = sin(pkin(12));
t24 = cos(pkin(12));
t16 = -t21 * t29 + t24 * t40;
t34 = t16 * t25 - t24 * t48;
t18 = -t21 * t40 - t24 * t29;
t33 = t18 * t25 + t21 * t48;
t19 = -t21 * t41 + t24 * t32;
t17 = t21 * t32 + t24 * t41;
t15 = t26 * t25 - t32 * t48;
t12 = -t18 * t22 + t21 * t44;
t11 = -t16 * t22 - t24 * t44;
t10 = t28 * t47 + (t25 * t38 + t37) * t23;
t8 = t18 * t31 - t19 * t43;
t6 = t16 * t31 - t17 * t43;
t4 = t19 * t31 + t33 * t28;
t2 = t17 * t31 + t34 * t28;
t1 = [0 (t19 * t46 + t8 * t30) * r_i_i_C(1) + (t19 * t45 - t8 * t27) * r_i_i_C(2) + t8 * pkin(3) + t18 * pkin(2) + t19 * t49 + t50 * (t18 * t28 + t19 * t42) t50 * t4 + t35 * (-t19 * t28 + t33 * t31) (t12 * t30 - t4 * t27) * r_i_i_C(1) + (-t12 * t27 - t4 * t30) * r_i_i_C(2), 0, 0; 0 (t17 * t46 + t6 * t30) * r_i_i_C(1) + (t17 * t45 - t6 * t27) * r_i_i_C(2) + t6 * pkin(3) + t16 * pkin(2) + t17 * t49 + t50 * (t16 * t28 + t17 * t42) t50 * t2 + t35 * (-t17 * t28 + t34 * t31) (t11 * t30 - t2 * t27) * r_i_i_C(1) + (-t11 * t27 - t2 * t30) * r_i_i_C(2), 0, 0; 1 (t35 * (-t25 * t39 + t36) + t32 * pkin(2) + (t27 * r_i_i_C(1) + t30 * r_i_i_C(2) + pkin(9)) * t29 * t22 + t50 * (t25 * t37 + t38)) * t23, t50 * t10 + t35 * (t31 * t47 + (t25 * t36 - t39) * t23) (-t10 * t27 + t15 * t30) * r_i_i_C(1) + (-t10 * t30 - t15 * t27) * r_i_i_C(2), 0, 0;];
Ja_transl  = t1;
