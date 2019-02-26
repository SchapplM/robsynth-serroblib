% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRP4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:47
% EndTime: 2019-02-26 21:26:47
% DurationCPUTime: 0.13s
% Computational Cost: add. (88->34), mult. (217->55), div. (0->0), fcn. (263->8), ass. (0->25)
t14 = cos(qJ(2));
t11 = sin(qJ(2));
t21 = -r_i_i_C(3) - pkin(8) + qJ(3);
t19 = t21 * t11;
t26 = t14 * pkin(2) + pkin(1) + t19;
t10 = sin(qJ(5));
t13 = cos(qJ(5));
t24 = pkin(3) + pkin(4);
t17 = t13 * r_i_i_C(1) - t10 * r_i_i_C(2) + t24;
t18 = t10 * r_i_i_C(1) + t13 * r_i_i_C(2) + qJ(4);
t8 = sin(pkin(9));
t9 = cos(pkin(9));
t25 = t17 * t9 + t18 * t8 + pkin(2);
t12 = sin(qJ(1));
t23 = t12 * t14;
t15 = cos(qJ(1));
t22 = t14 * t15;
t16 = -t25 * t11 + t21 * t14;
t6 = t12 * t8 + t9 * t22;
t5 = -t12 * t9 + t8 * t22;
t4 = -t15 * t8 + t9 * t23;
t3 = t15 * t9 + t8 * t23;
t2 = t5 * t10 + t6 * t13;
t1 = -t6 * t10 + t5 * t13;
t7 = [t15 * pkin(7) - t26 * t12 - t17 * t4 - t18 * t3, t16 * t15, t15 * t11, t5, r_i_i_C(1) * t1 - t2 * r_i_i_C(2), 0; t12 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t5 * qJ(4) + t26 * t15 + t24 * t6, t16 * t12, t12 * t11, t3 (-t4 * t10 + t3 * t13) * r_i_i_C(1) + (-t3 * t10 - t4 * t13) * r_i_i_C(2), 0; 0, t25 * t14 + t19, -t14, t11 * t8 ((-t10 * t9 + t13 * t8) * r_i_i_C(1) + (-t10 * t8 - t13 * t9) * r_i_i_C(2)) * t11, 0;];
Ja_transl  = t7;
