% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP8_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP8_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:13:08
% EndTime: 2019-02-26 22:13:09
% DurationCPUTime: 0.15s
% Computational Cost: add. (119->39), mult. (291->61), div. (0->0), fcn. (354->8), ass. (0->30)
t18 = cos(qJ(2));
t14 = sin(qJ(2));
t28 = -r_i_i_C(3) - pkin(9) + pkin(8);
t25 = t28 * t14;
t34 = t18 * pkin(2) + pkin(1) + t25;
t13 = sin(qJ(3));
t17 = cos(qJ(3));
t16 = cos(qJ(5));
t12 = sin(qJ(5));
t32 = pkin(3) + pkin(4);
t24 = t12 * r_i_i_C(2) - t32;
t21 = t16 * r_i_i_C(1) - t24;
t26 = -t12 * r_i_i_C(1) - qJ(4);
t22 = t16 * r_i_i_C(2) - t26;
t33 = t22 * t13 + t21 * t17 + pkin(2);
t19 = cos(qJ(1));
t29 = t19 * t17;
t15 = sin(qJ(1));
t31 = t15 * t18;
t6 = t13 * t31 + t29;
t3 = t6 * t16;
t30 = t19 * t13;
t23 = (-(t12 * t17 - t13 * t16) * r_i_i_C(1) - (t12 * t13 + t16 * t17) * r_i_i_C(2)) * t14;
t20 = -t33 * t14 + t28 * t18;
t9 = t15 * t13 + t18 * t29;
t8 = -t15 * t17 + t18 * t30;
t7 = t17 * t31 - t30;
t2 = t8 * t12 + t9 * t16;
t1 = -t9 * t12 + t8 * t16;
t4 = [t19 * pkin(7) - t3 * r_i_i_C(2) - t34 * t15 - t21 * t7 + t26 * t6, t20 * t19, -t21 * t8 + t22 * t9, t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t15 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t8 * qJ(4) + t34 * t19 + t32 * t9, t20 * t15, -t3 * r_i_i_C(1) + t22 * t7 + t24 * t6, t6 (-t7 * t12 + t3) * r_i_i_C(1) + (-t6 * t12 - t7 * t16) * r_i_i_C(2), 0; 0, t33 * t18 + t25 (qJ(4) * t17 - t32 * t13) * t14 - t23, t14 * t13, t23, 0;];
Ja_transl  = t4;
