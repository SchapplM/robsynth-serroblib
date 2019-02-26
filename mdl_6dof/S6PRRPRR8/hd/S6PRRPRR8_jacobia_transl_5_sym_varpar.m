% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR8_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR8_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:07
% EndTime: 2019-02-26 20:08:07
% DurationCPUTime: 0.16s
% Computational Cost: add. (206->52), mult. (575->98), div. (0->0), fcn. (742->12), ass. (0->40)
t26 = sin(pkin(7));
t27 = sin(pkin(6));
t54 = t26 * t27;
t29 = cos(pkin(7));
t53 = t27 * t29;
t31 = sin(qJ(3));
t52 = t29 * t31;
t34 = cos(qJ(3));
t51 = t29 * t34;
t32 = sin(qJ(2));
t50 = t31 * t32;
t35 = cos(qJ(2));
t49 = t31 * t35;
t48 = t32 * t34;
t47 = t34 * t35;
t46 = cos(pkin(6));
t45 = sin(pkin(12));
t44 = r_i_i_C(3) + pkin(10) + pkin(3);
t28 = cos(pkin(12));
t43 = t28 * t54;
t42 = t27 * t45;
t41 = t28 * t46;
t40 = t46 * t26;
t39 = t26 * t42;
t38 = t46 * t45;
t30 = sin(qJ(5));
t33 = cos(qJ(5));
t37 = t30 * r_i_i_C(1) + t33 * r_i_i_C(2) + qJ(4);
t36 = (t33 * r_i_i_C(1) - t30 * r_i_i_C(2) + pkin(4) + pkin(9)) * t26;
t21 = t28 * t35 - t32 * t38;
t20 = -t28 * t32 - t35 * t38;
t19 = t32 * t41 + t45 * t35;
t18 = -t45 * t32 + t35 * t41;
t17 = t46 * t29 - t35 * t54;
t12 = -t20 * t26 + t29 * t42;
t11 = -t18 * t26 - t28 * t53;
t9 = t27 * t50 - t34 * t40 - t47 * t53;
t3 = -t20 * t51 + t21 * t31 - t34 * t39;
t1 = -t18 * t51 + t19 * t31 + t34 * t43;
t2 = [0, t20 * pkin(2) + t37 * (t20 * t31 + t21 * t51) + t21 * t36 + t44 * (t20 * t34 - t21 * t52) t37 * (t21 * t34 + (t20 * t29 + t39) * t31) - t44 * t3, t3 (-t12 * t30 + t3 * t33) * r_i_i_C(1) + (-t12 * t33 - t3 * t30) * r_i_i_C(2), 0; 0, t18 * pkin(2) + t37 * (t18 * t31 + t19 * t51) + t19 * t36 + t44 * (t18 * t34 - t19 * t52) t37 * (t19 * t34 + (t18 * t29 - t43) * t31) - t44 * t1, t1 (t1 * t33 - t11 * t30) * r_i_i_C(1) + (-t1 * t30 - t11 * t33) * r_i_i_C(2), 0; 1 (t37 * (t29 * t48 + t49) + t35 * pkin(2) + t32 * t36 + t44 * (-t29 * t50 + t47)) * t27, -t44 * t9 + t37 * (t31 * t40 + (t29 * t49 + t48) * t27) t9 (-t17 * t30 + t9 * t33) * r_i_i_C(1) + (-t17 * t33 - t9 * t30) * r_i_i_C(2), 0;];
Ja_transl  = t2;
