% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:37:35
% EndTime: 2019-02-26 20:37:35
% DurationCPUTime: 0.11s
% Computational Cost: add. (98->30), mult. (128->43), div. (0->0), fcn. (144->8), ass. (0->27)
t14 = sin(qJ(4));
t17 = cos(qJ(4));
t28 = r_i_i_C(3) + pkin(9) + pkin(8);
t25 = t28 * t17;
t16 = cos(qJ(5));
t9 = pkin(5) * t16 + pkin(4);
t32 = -t14 * t9 - pkin(1) - qJ(3) + t25;
t12 = qJ(5) + qJ(6);
t10 = sin(t12);
t11 = cos(t12);
t18 = cos(qJ(1));
t15 = sin(qJ(1));
t27 = t14 * t15;
t5 = t10 * t27 - t11 * t18;
t6 = -t10 * t18 - t11 * t27;
t31 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t26 = t14 * t18;
t7 = -t10 * t26 - t11 * t15;
t8 = -t10 * t15 + t11 * t26;
t30 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t13 = sin(qJ(5));
t29 = pkin(5) * t13;
t23 = -pkin(7) + qJ(2) - t29;
t22 = -r_i_i_C(1) * t10 - r_i_i_C(2) * t11;
t21 = r_i_i_C(1) * t11 - r_i_i_C(2) * t10 + t9;
t20 = t28 * t14 + t21 * t17;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t32 * t15 + t23 * t18, t15, t18, t20 * t18 (-t13 * t26 - t15 * t16) * pkin(5) + t30, t30; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t23 * t15 - t32 * t18, -t18, t15, t20 * t15 (-t13 * t27 + t16 * t18) * pkin(5) + t31, t31; 0, 0, 0, -t21 * t14 + t25 (t22 - t29) * t17, t22 * t17;];
Ja_transl  = t1;
