% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR7_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR7_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:47
% EndTime: 2019-02-26 19:49:47
% DurationCPUTime: 0.09s
% Computational Cost: add. (74->23), mult. (186->40), div. (0->0), fcn. (235->8), ass. (0->23)
t28 = pkin(4) - r_i_i_C(2);
t13 = sin(pkin(6));
t16 = sin(qJ(4));
t27 = t13 * t16;
t18 = cos(qJ(4));
t26 = t13 * t18;
t19 = cos(qJ(2));
t25 = t13 * t19;
t15 = cos(pkin(6));
t17 = sin(qJ(2));
t24 = t15 * t17;
t23 = t15 * t19;
t22 = r_i_i_C(3) + qJ(5);
t21 = pkin(2) + pkin(8) + r_i_i_C(1);
t20 = t28 * t16 - t22 * t18 + qJ(3);
t14 = cos(pkin(10));
t12 = sin(pkin(10));
t9 = t15 * t16 + t18 * t25;
t7 = t12 * t23 + t14 * t17;
t5 = t12 * t17 - t14 * t23;
t3 = t14 * t27 + t18 * t5;
t1 = t12 * t27 - t18 * t7;
t2 = [0, -t21 * t7 + t20 * (-t12 * t24 + t14 * t19) t7, t22 * (t12 * t26 + t16 * t7) - t28 * t1, t1, 0; 0, -t21 * t5 + t20 * (t12 * t19 + t14 * t24) t5, -t22 * (t14 * t26 - t16 * t5) + t28 * t3, -t3, 0; 1 (t20 * t17 + t21 * t19) * t13, -t25, -t28 * t9 + t22 * (t15 * t18 - t16 * t25) t9, 0;];
Ja_transl  = t2;
