% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:38
% EndTime: 2019-02-26 19:47:38
% DurationCPUTime: 0.18s
% Computational Cost: add. (241->41), mult. (625->75), div. (0->0), fcn. (824->12), ass. (0->33)
t24 = sin(pkin(11));
t31 = sin(qJ(2));
t34 = cos(qJ(2));
t43 = cos(pkin(11));
t37 = -t31 * t24 + t34 * t43;
t30 = sin(qJ(4));
t33 = cos(qJ(4));
t29 = sin(qJ(6));
t32 = cos(qJ(6));
t39 = t29 * r_i_i_C(1) + t32 * r_i_i_C(2) + qJ(5);
t42 = pkin(4) + pkin(9) + r_i_i_C(3);
t35 = t39 * t30 + t42 * t33 + pkin(3);
t25 = sin(pkin(10));
t26 = sin(pkin(6));
t47 = t25 * t26;
t27 = cos(pkin(10));
t46 = t27 * t26;
t28 = cos(pkin(6));
t45 = t28 * t34;
t19 = -t34 * t24 - t31 * t43;
t17 = t19 * t28;
t7 = -t27 * t17 + t25 * t37;
t40 = -t25 * t17 - t27 * t37;
t38 = t32 * r_i_i_C(1) - t29 * r_i_i_C(2) + pkin(5) + pkin(8);
t36 = t37 * t28;
t16 = t19 * t26;
t15 = t37 * t26;
t11 = -t16 * t30 - t28 * t33;
t9 = t27 * t19 - t25 * t36;
t6 = t25 * t19 + t27 * t36;
t3 = -t30 * t40 - t33 * t47;
t1 = t7 * t30 + t33 * t46;
t2 = [0 (-t25 * t45 - t27 * t31) * pkin(2) - t38 * t40 + t35 * t9, t47, t39 * (t30 * t47 - t33 * t40) - t42 * t3, t3 (t9 * t29 + t3 * t32) * r_i_i_C(1) + (-t3 * t29 + t9 * t32) * r_i_i_C(2); 0 (-t25 * t31 + t27 * t45) * pkin(2) + t38 * t7 + t35 * t6, -t46, t39 * (-t30 * t46 + t7 * t33) - t42 * t1, t1 (t1 * t32 + t6 * t29) * r_i_i_C(1) + (-t1 * t29 + t6 * t32) * r_i_i_C(2); 1, t26 * t34 * pkin(2) + t35 * t15 - t38 * t16, t28, t39 * (-t16 * t33 + t28 * t30) - t42 * t11, t11 (t11 * t32 + t15 * t29) * r_i_i_C(1) + (-t11 * t29 + t15 * t32) * r_i_i_C(2);];
Ja_transl  = t2;
