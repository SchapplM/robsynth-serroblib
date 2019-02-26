% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPPRR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:23:46
% EndTime: 2019-02-26 20:23:46
% DurationCPUTime: 0.11s
% Computational Cost: add. (115->32), mult. (167->41), div. (0->0), fcn. (217->9), ass. (0->25)
t11 = pkin(10) + qJ(5);
t10 = cos(t11);
t13 = sin(qJ(6));
t14 = cos(qJ(6));
t16 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + pkin(5);
t25 = pkin(8) + r_i_i_C(3);
t9 = sin(t11);
t18 = t25 * t9;
t28 = -t16 * t10 - t18;
t27 = pkin(1) + pkin(2);
t24 = cos(qJ(1));
t23 = sin(qJ(1));
t22 = t10 * t13;
t21 = t10 * t14;
t20 = cos(pkin(9));
t19 = sin(pkin(9));
t17 = t13 * r_i_i_C(1) + t14 * r_i_i_C(2);
t15 = -t25 * t10 + t16 * t9;
t12 = -pkin(7) - qJ(4);
t8 = cos(pkin(10)) * pkin(4) + pkin(3);
t4 = t24 * t19 - t23 * t20;
t3 = -t23 * t19 - t24 * t20;
t2 = t4 * t13 - t3 * t21;
t1 = t14 * t4 + t3 * t22;
t5 = [t24 * qJ(2) + (-t12 + t17) * t3 + (t8 - t28) * t4 - t27 * t23, t23, 0, t4, t15 * t3, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t23 * qJ(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t4 * t12 + (-t10 * pkin(5) - t18 - t8) * t3 + t27 * t24, -t24, 0, -t3, t15 * t4 (-t14 * t3 + t4 * t22) * r_i_i_C(1) + (t3 * t13 + t4 * t21) * r_i_i_C(2); 0, 0, -1, 0, t28, t17 * t9;];
Ja_transl  = t5;
