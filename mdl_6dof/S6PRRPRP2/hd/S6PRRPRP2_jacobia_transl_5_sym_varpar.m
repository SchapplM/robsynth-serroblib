% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRP2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:47
% EndTime: 2019-02-26 20:01:47
% DurationCPUTime: 0.16s
% Computational Cost: add. (176->42), mult. (300->74), div. (0->0), fcn. (376->12), ass. (0->32)
t15 = qJ(3) + pkin(11);
t13 = sin(t15);
t14 = cos(t15);
t25 = cos(qJ(3));
t21 = sin(qJ(5));
t24 = cos(qJ(5));
t29 = t24 * r_i_i_C(1) - t21 * r_i_i_C(2) + pkin(4);
t37 = pkin(9) + r_i_i_C(3);
t38 = t25 * pkin(3) + t37 * t13 + t29 * t14 + pkin(2);
t16 = sin(pkin(10));
t17 = sin(pkin(6));
t36 = t16 * t17;
t18 = cos(pkin(10));
t35 = t17 * t18;
t23 = sin(qJ(2));
t34 = t17 * t23;
t33 = t17 * t25;
t26 = cos(qJ(2));
t32 = t17 * t26;
t19 = cos(pkin(6));
t31 = t19 * t23;
t30 = t19 * t26;
t28 = t21 * r_i_i_C(1) + t24 * r_i_i_C(2) + pkin(8) + qJ(4);
t22 = sin(qJ(3));
t10 = -t16 * t31 + t18 * t26;
t9 = t16 * t30 + t18 * t23;
t8 = t16 * t26 + t18 * t31;
t7 = t16 * t23 - t18 * t30;
t6 = t19 * t13 + t14 * t34;
t4 = t10 * t14 + t13 * t36;
t2 = -t13 * t35 + t8 * t14;
t1 = [0, t28 * t10 - t38 * t9, t37 * t4 + (-t10 * t22 + t16 * t33) * pkin(3) + t29 * (-t10 * t13 + t14 * t36) t9 (-t4 * t21 + t9 * t24) * r_i_i_C(1) + (-t9 * t21 - t4 * t24) * r_i_i_C(2), 0; 0, t28 * t8 - t38 * t7, t37 * t2 + (-t18 * t33 - t22 * t8) * pkin(3) + t29 * (-t8 * t13 - t14 * t35) t7 (-t2 * t21 + t7 * t24) * r_i_i_C(1) + (-t2 * t24 - t7 * t21) * r_i_i_C(2), 0; 1 (t28 * t23 + t38 * t26) * t17, t37 * t6 + (t19 * t25 - t22 * t34) * pkin(3) + t29 * (-t13 * t34 + t19 * t14) -t32 (-t6 * t21 - t24 * t32) * r_i_i_C(1) + (t21 * t32 - t6 * t24) * r_i_i_C(2), 0;];
Ja_transl  = t1;
