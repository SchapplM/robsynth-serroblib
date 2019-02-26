% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR3
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
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR3_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR3_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobia_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:23
% EndTime: 2019-02-26 20:05:24
% DurationCPUTime: 0.15s
% Computational Cost: add. (112->44), mult. (301->84), div. (0->0), fcn. (383->12), ass. (0->32)
t25 = cos(qJ(3));
t36 = t25 * pkin(3);
t16 = sin(pkin(12));
t18 = sin(pkin(6));
t35 = t16 * t18;
t17 = sin(pkin(7));
t34 = t17 * t18;
t20 = cos(pkin(12));
t33 = t18 * t20;
t21 = cos(pkin(7));
t32 = t18 * t21;
t22 = cos(pkin(6));
t24 = sin(qJ(2));
t31 = t22 * t24;
t26 = cos(qJ(2));
t30 = t22 * t26;
t15 = sin(pkin(13));
t19 = cos(pkin(13));
t23 = sin(qJ(3));
t29 = t25 * t15 + t23 * t19;
t10 = t23 * t15 - t25 * t19;
t28 = -t10 * r_i_i_C(1) - r_i_i_C(2) * t29 + pkin(2) + t36;
t3 = t10 * t21;
t4 = t29 * t21;
t27 = -t21 * t23 * pkin(3) - t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + (r_i_i_C(3) + pkin(9) + qJ(4)) * t17;
t9 = -t16 * t31 + t20 * t26;
t8 = -t16 * t30 - t20 * t24;
t7 = t16 * t26 + t20 * t31;
t6 = -t16 * t24 + t20 * t30;
t2 = t29 * t17;
t1 = t10 * t17;
t5 = [0, t27 * t9 + t28 * t8 (-t1 * t35 - t29 * t9 - t8 * t3) * r_i_i_C(1) + (t9 * t10 - t2 * t35 - t8 * t4) * r_i_i_C(2) + (-t9 * t23 + (t16 * t34 + t21 * t8) * t25) * pkin(3), t16 * t32 - t8 * t17, 0, 0; 0, t27 * t7 + t28 * t6 (t1 * t33 - t29 * t7 - t6 * t3) * r_i_i_C(1) + (t7 * t10 + t2 * t33 - t6 * t4) * r_i_i_C(2) + (-t7 * t23 + (-t17 * t33 + t21 * t6) * t25) * pkin(3), -t6 * t17 - t20 * t32, 0, 0; 1 (t27 * t24 + t28 * t26) * t18 (-t1 * r_i_i_C(1) - t2 * r_i_i_C(2) + t17 * t36) * t22 + ((-t24 * t29 - t26 * t3) * r_i_i_C(1) + (t10 * t24 - t26 * t4) * r_i_i_C(2) + (t26 * t21 * t25 - t24 * t23) * pkin(3)) * t18, t22 * t21 - t26 * t34, 0, 0;];
Ja_transl  = t5;
