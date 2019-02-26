% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR10_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR10_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:58
% EndTime: 2019-02-26 20:53:59
% DurationCPUTime: 0.11s
% Computational Cost: add. (155->35), mult. (135->47), div. (0->0), fcn. (152->10), ass. (0->29)
t20 = cos(qJ(3));
t29 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(4);
t34 = t29 * t20;
t18 = sin(qJ(3));
t17 = pkin(10) + qJ(5);
t15 = qJ(6) + t17;
t11 = sin(t15);
t12 = cos(t15);
t14 = cos(t17);
t9 = pkin(5) * t14 + cos(pkin(10)) * pkin(4) + pkin(3);
t23 = r_i_i_C(1) * t12 - r_i_i_C(2) * t11 + t9;
t33 = t29 * t18 + t23 * t20;
t21 = cos(qJ(1));
t19 = sin(qJ(1));
t27 = t19 * t11;
t5 = t12 * t21 - t18 * t27;
t26 = t19 * t12;
t6 = t11 * t21 + t18 * t26;
t32 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t28 = t18 * t21;
t7 = t11 * t28 + t26;
t8 = t12 * t28 - t27;
t31 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
t13 = sin(t17);
t30 = pkin(5) * t13;
t25 = pkin(1) + pkin(7) + t30 + sin(pkin(10)) * pkin(4);
t24 = -r_i_i_C(1) * t11 - r_i_i_C(2) * t12;
t22 = t18 * t9 + qJ(2) - t34;
t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) - t25 * t19 + t22 * t21, t19, t33 * t19, -t19 * t20 (-t13 * t18 * t19 + t14 * t21) * pkin(5) + t32, t32; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t22 * t19 + t25 * t21, -t21, -t33 * t21, t21 * t20 (t13 * t28 + t14 * t19) * pkin(5) + t31, t31; 0, 0, -t23 * t18 + t34, t18 (t24 - t30) * t20, t24 * t20;];
Ja_transl  = t1;
