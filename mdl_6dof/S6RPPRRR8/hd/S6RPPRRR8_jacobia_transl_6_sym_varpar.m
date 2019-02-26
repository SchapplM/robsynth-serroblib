% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRR8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:38:39
% EndTime: 2019-02-26 20:38:39
% DurationCPUTime: 0.10s
% Computational Cost: add. (140->32), mult. (130->43), div. (0->0), fcn. (146->9), ass. (0->31)
t14 = pkin(10) + qJ(4);
t11 = cos(t14);
t32 = r_i_i_C(3) + pkin(9) + pkin(8);
t37 = t32 * t11;
t10 = sin(t14);
t15 = qJ(5) + qJ(6);
t12 = sin(t15);
t13 = cos(t15);
t20 = cos(qJ(5));
t9 = pkin(5) * t20 + pkin(4);
t24 = r_i_i_C(1) * t13 - r_i_i_C(2) * t12 + t9;
t36 = t32 * t10 + t24 * t11;
t21 = cos(qJ(1));
t27 = t13 * t21;
t19 = sin(qJ(1));
t30 = t12 * t19;
t5 = -t10 * t30 + t27;
t28 = t13 * t19;
t29 = t12 * t21;
t6 = t10 * t28 + t29;
t35 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t7 = t10 * t29 + t28;
t8 = t10 * t27 - t30;
t34 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
t18 = sin(qJ(5));
t33 = pkin(5) * t18;
t31 = t10 * t18;
t26 = pkin(1) + pkin(7) + qJ(3) + t33;
t25 = -r_i_i_C(1) * t12 - r_i_i_C(2) * t13;
t23 = pkin(3) * sin(pkin(10)) + t10 * t9 - t37 + qJ(2);
t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) - t26 * t19 + t23 * t21, t19, t21, t36 * t19 (-t19 * t31 + t20 * t21) * pkin(5) + t35, t35; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t23 * t19 + t26 * t21, -t21, t19, -t36 * t21 (t19 * t20 + t21 * t31) * pkin(5) + t34, t34; 0, 0, 0, -t24 * t10 + t37 (t25 - t33) * t11, t25 * t11;];
Ja_transl  = t1;
