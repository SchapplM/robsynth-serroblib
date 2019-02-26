% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:38
% EndTime: 2019-02-26 20:10:38
% DurationCPUTime: 0.11s
% Computational Cost: add. (160->37), mult. (189->60), div. (0->0), fcn. (228->12), ass. (0->29)
t23 = qJ(3) + qJ(4);
t18 = pkin(12) + t23;
t15 = sin(t18);
t16 = cos(t18);
t25 = sin(pkin(6));
t26 = cos(pkin(11));
t34 = t25 * t26;
t24 = sin(pkin(11));
t29 = cos(qJ(2));
t27 = cos(pkin(6));
t28 = sin(qJ(2));
t32 = t27 * t28;
t8 = t24 * t29 + t26 * t32;
t39 = (-t8 * t15 - t16 * t34) * r_i_i_C(1) + (t15 * t34 - t8 * t16) * r_i_i_C(2);
t10 = -t24 * t32 + t26 * t29;
t35 = t24 * t25;
t38 = (-t10 * t15 + t16 * t35) * r_i_i_C(1) + (-t10 * t16 - t15 * t35) * r_i_i_C(2);
t33 = t25 * t28;
t37 = (-t15 * t33 + t27 * t16) * r_i_i_C(1) + (-t27 * t15 - t16 * t33) * r_i_i_C(2);
t36 = r_i_i_C(3) + qJ(5) + pkin(9) + pkin(8);
t31 = t27 * t29;
t20 = cos(t23);
t13 = pkin(4) * t20 + cos(qJ(3)) * pkin(3);
t30 = r_i_i_C(1) * t16 - r_i_i_C(2) * t15 + pkin(2) + t13;
t19 = sin(t23);
t12 = -pkin(4) * t19 - sin(qJ(3)) * pkin(3);
t9 = t24 * t31 + t26 * t28;
t7 = t24 * t28 - t26 * t31;
t1 = [0, t36 * t10 - t30 * t9, t10 * t12 + t13 * t35 + t38 (-t10 * t19 + t20 * t35) * pkin(4) + t38, t9, 0; 0, -t30 * t7 + t36 * t8, t8 * t12 - t13 * t34 + t39 (-t19 * t8 - t20 * t34) * pkin(4) + t39, t7, 0; 1 (t36 * t28 + t30 * t29) * t25, t12 * t33 + t27 * t13 + t37 (-t19 * t33 + t20 * t27) * pkin(4) + t37, -t25 * t29, 0;];
Ja_transl  = t1;
