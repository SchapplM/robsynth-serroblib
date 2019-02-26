% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPPRR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:24:53
% EndTime: 2019-02-26 20:24:53
% DurationCPUTime: 0.10s
% Computational Cost: add. (80->30), mult. (163->40), div. (0->0), fcn. (211->8), ass. (0->23)
t23 = pkin(1) + qJ(3);
t22 = pkin(3) + qJ(2);
t10 = cos(qJ(5));
t7 = sin(qJ(6));
t9 = cos(qJ(6));
t13 = r_i_i_C(1) * t9 - r_i_i_C(2) * t7 + pkin(5);
t20 = pkin(8) + r_i_i_C(3);
t8 = sin(qJ(5));
t21 = t20 * t10 - t13 * t8;
t19 = t10 * t7;
t18 = t10 * t9;
t17 = sin(qJ(1));
t16 = sin(pkin(9));
t15 = t20 * t8;
t14 = -t7 * r_i_i_C(1) - t9 * r_i_i_C(2);
t12 = t13 * t10 + t15;
t11 = cos(qJ(1));
t6 = cos(pkin(9));
t4 = t11 * t16 + t17 * t6;
t3 = t11 * t6 - t17 * t16;
t2 = t4 * t18 - t3 * t7;
t1 = -t4 * t19 - t3 * t9;
t5 = [t22 * t11 + (pkin(7) - t14) * t4 + (pkin(4) + t12) * t3 - t23 * t17, t17, t11, 0, t21 * t4, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; -t3 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t23 * t11 + (pkin(5) * t10 + pkin(4) + t15) * t4 + t22 * t17, -t11, t17, 0, -t21 * t3 (t3 * t19 - t4 * t9) * r_i_i_C(1) + (t3 * t18 + t4 * t7) * r_i_i_C(2); 0, 0, 0, 1, t12, t14 * t8;];
Ja_transl  = t5;
