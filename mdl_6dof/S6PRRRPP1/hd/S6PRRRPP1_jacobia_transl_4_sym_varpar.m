% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPP1_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPP1_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:54
% EndTime: 2019-02-26 20:08:54
% DurationCPUTime: 0.12s
% Computational Cost: add. (102->32), mult. (269->63), div. (0->0), fcn. (338->10), ass. (0->28)
t15 = sin(qJ(3));
t18 = cos(qJ(3));
t14 = sin(qJ(4));
t17 = cos(qJ(4));
t22 = t17 * r_i_i_C(1) - t14 * r_i_i_C(2) + pkin(3);
t31 = pkin(9) + r_i_i_C(3);
t32 = t31 * t15 + t22 * t18 + pkin(2);
t13 = sin(pkin(6));
t30 = t13 * t15;
t29 = t13 * t18;
t19 = cos(qJ(2));
t28 = t13 * t19;
t27 = cos(pkin(6));
t26 = cos(pkin(10));
t12 = sin(pkin(10));
t25 = t12 * t27;
t24 = t13 * t26;
t23 = t27 * t26;
t21 = t14 * r_i_i_C(1) + t17 * r_i_i_C(2) + pkin(8);
t16 = sin(qJ(2));
t10 = t27 * t15 + t16 * t29;
t8 = -t16 * t25 + t26 * t19;
t7 = t26 * t16 + t19 * t25;
t6 = t12 * t19 + t16 * t23;
t5 = t12 * t16 - t19 * t23;
t4 = t12 * t30 + t8 * t18;
t2 = -t15 * t24 + t6 * t18;
t1 = [0, t21 * t8 - t32 * t7, t31 * t4 + t22 * (t12 * t29 - t8 * t15) (-t4 * t14 + t7 * t17) * r_i_i_C(1) + (-t7 * t14 - t4 * t17) * r_i_i_C(2), 0, 0; 0, t21 * t6 - t32 * t5, t31 * t2 + t22 * (-t6 * t15 - t18 * t24) (-t2 * t14 + t5 * t17) * r_i_i_C(1) + (-t5 * t14 - t2 * t17) * r_i_i_C(2), 0, 0; 1 (t21 * t16 + t32 * t19) * t13, t31 * t10 + t22 * (-t16 * t30 + t27 * t18) (-t10 * t14 - t17 * t28) * r_i_i_C(1) + (-t10 * t17 + t14 * t28) * r_i_i_C(2), 0, 0;];
Ja_transl  = t1;
