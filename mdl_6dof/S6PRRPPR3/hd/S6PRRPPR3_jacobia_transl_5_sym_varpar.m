% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPPR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:59:23
% EndTime: 2019-02-26 19:59:23
% DurationCPUTime: 0.09s
% Computational Cost: add. (86->21), mult. (217->40), div. (0->0), fcn. (273->8), ass. (0->24)
t15 = sin(pkin(6));
t18 = sin(qJ(3));
t29 = t15 * t18;
t20 = cos(qJ(3));
t28 = t15 * t20;
t17 = cos(pkin(6));
t19 = sin(qJ(2));
t27 = t17 * t19;
t21 = cos(qJ(2));
t26 = t17 * t21;
t25 = r_i_i_C(1) + qJ(4);
t24 = pkin(3) + pkin(4) - r_i_i_C(2);
t23 = pkin(8) - r_i_i_C(3) - qJ(5);
t22 = t25 * t18 + t24 * t20 + pkin(2);
t16 = cos(pkin(10));
t14 = sin(pkin(10));
t9 = -t17 * t20 + t19 * t29;
t8 = -t14 * t27 + t16 * t21;
t7 = -t14 * t26 - t16 * t19;
t6 = t14 * t21 + t16 * t27;
t5 = -t14 * t19 + t16 * t26;
t3 = -t14 * t28 + t18 * t8;
t1 = t16 * t28 + t18 * t6;
t2 = [0, t22 * t7 + t23 * t8, t25 * (t14 * t29 + t20 * t8) - t24 * t3, t3, t7, 0; 0, t22 * t5 + t23 * t6, t25 * (-t16 * t29 + t20 * t6) - t24 * t1, t1, t5, 0; 1 (t23 * t19 + t22 * t21) * t15, t25 * (t17 * t18 + t19 * t28) - t24 * t9, t9, t15 * t21, 0;];
Ja_transl  = t2;
