% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR6
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR6_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR6_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobia_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:05
% EndTime: 2019-02-26 20:07:05
% DurationCPUTime: 0.14s
% Computational Cost: add. (139->47), mult. (392->87), div. (0->0), fcn. (504->12), ass. (0->39)
t24 = sin(pkin(7));
t51 = t24 * pkin(9);
t23 = sin(pkin(13));
t50 = t23 * t24;
t25 = sin(pkin(6));
t49 = t24 * t25;
t26 = cos(pkin(13));
t48 = t24 * t26;
t28 = cos(pkin(7));
t29 = sin(qJ(3));
t47 = t28 * t29;
t31 = cos(qJ(3));
t46 = t28 * t31;
t30 = sin(qJ(2));
t45 = t29 * t30;
t32 = cos(qJ(2));
t44 = t29 * t32;
t43 = t30 * t31;
t42 = t31 * t32;
t41 = r_i_i_C(3) + qJ(4);
t40 = cos(pkin(6));
t39 = sin(pkin(12));
t27 = cos(pkin(12));
t38 = t27 * t49;
t37 = t27 * t40;
t36 = t40 * t24;
t35 = t39 * t49;
t34 = t40 * t39;
t33 = t26 * r_i_i_C(1) - t23 * r_i_i_C(2) + pkin(3);
t18 = t27 * t32 - t30 * t34;
t17 = -t27 * t30 - t32 * t34;
t16 = t30 * t37 + t39 * t32;
t15 = -t39 * t30 + t32 * t37;
t9 = -t31 * t36 + (-t28 * t42 + t45) * t25;
t8 = t17 * t31 - t18 * t47;
t6 = t15 * t31 - t16 * t47;
t3 = -t17 * t46 + t18 * t29 - t31 * t35;
t1 = -t15 * t46 + t16 * t29 + t31 * t38;
t2 = [0 (t18 * t50 + t8 * t26) * r_i_i_C(1) + (t18 * t48 - t8 * t23) * r_i_i_C(2) + t8 * pkin(3) + t17 * pkin(2) + t18 * t51 + t41 * (t17 * t29 + t18 * t46) t41 * (t18 * t31 + (t17 * t28 + t35) * t29) - t33 * t3, t3, 0, 0; 0 (t16 * t50 + t6 * t26) * r_i_i_C(1) + (t16 * t48 - t6 * t23) * r_i_i_C(2) + t6 * pkin(3) + t15 * pkin(2) + t16 * t51 + t41 * (t15 * t29 + t16 * t46) t41 * (t16 * t31 + (t15 * t28 - t38) * t29) - t33 * t1, t1, 0, 0; 1 (t33 * (-t28 * t45 + t42) + t32 * pkin(2) + (r_i_i_C(1) * t23 + r_i_i_C(2) * t26 + pkin(9)) * t30 * t24 + t41 * (t28 * t43 + t44)) * t25, t41 * (t29 * t36 + (t28 * t44 + t43) * t25) - t33 * t9, t9, 0, 0;];
Ja_transl  = t2;
