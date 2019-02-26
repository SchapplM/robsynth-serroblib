% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPP1_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPP1_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobia_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:02:46
% EndTime: 2019-02-26 22:02:47
% DurationCPUTime: 0.17s
% Computational Cost: add. (115->36), mult. (318->59), div. (0->0), fcn. (385->10), ass. (0->32)
t10 = sin(qJ(3));
t13 = cos(qJ(3));
t6 = sin(pkin(10));
t8 = cos(pkin(10));
t22 = r_i_i_C(1) * t6 + r_i_i_C(2) * t8;
t24 = r_i_i_C(3) + qJ(4);
t7 = sin(pkin(6));
t9 = cos(pkin(6));
t18 = t22 * t9 - t24 * t7;
t21 = t8 * r_i_i_C(1) - t6 * r_i_i_C(2) + pkin(3);
t34 = t18 * t10 - t21 * t13 - pkin(2);
t11 = sin(qJ(2));
t20 = t22 * t7 + pkin(9);
t17 = t24 * t9 + t20;
t33 = t17 * t11;
t12 = sin(qJ(1));
t14 = cos(qJ(2));
t15 = cos(qJ(1));
t26 = t15 * t10;
t4 = t12 * t13 - t14 * t26;
t30 = t4 * t7;
t29 = t4 * t9;
t28 = t11 * t9;
t27 = t12 * t14;
t25 = t15 * t13;
t23 = -t14 * pkin(2) - pkin(1);
t16 = t34 * t11 + t17 * t14;
t5 = t12 * t10 + t14 * t25;
t3 = -t13 * t27 + t26;
t2 = t10 * t27 + t25;
t1 = t15 * t28 - t30;
t19 = [t15 * pkin(8) + t21 * t3 + t18 * t2 + (t23 - t33) * t12, t16 * t15, -t18 * t5 + t21 * t4, t1, 0, 0; (t6 * t29 + t5 * t8) * r_i_i_C(1) + (t8 * t29 - t5 * t6) * r_i_i_C(2) + t1 * r_i_i_C(3) + t5 * pkin(3) - qJ(4) * t30 + t12 * pkin(8) + ((t9 * qJ(4) + t20) * t11 - t23) * t15, t16 * t12, t18 * t3 - t21 * t2, t12 * t28 + t2 * t7, 0, 0; 0, -t34 * t14 + t33 (-t21 * t10 - t13 * t18) * t11, t11 * t10 * t7 - t14 * t9, 0, 0;];
Ja_transl  = t19;
