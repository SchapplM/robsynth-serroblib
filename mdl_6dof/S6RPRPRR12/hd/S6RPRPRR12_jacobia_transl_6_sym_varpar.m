% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR12_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR12_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:02
% EndTime: 2019-02-26 20:55:02
% DurationCPUTime: 0.11s
% Computational Cost: add. (104->32), mult. (146->43), div. (0->0), fcn. (163->8), ass. (0->28)
t15 = sin(qJ(3));
t24 = pkin(3) + r_i_i_C(3) + pkin(9) + pkin(8);
t33 = t24 * t15;
t18 = cos(qJ(3));
t13 = qJ(5) + qJ(6);
t11 = sin(t13);
t12 = cos(t13);
t14 = sin(qJ(5));
t23 = -pkin(5) * t14 - qJ(4);
t22 = r_i_i_C(1) * t11 + r_i_i_C(2) * t12 - t23;
t32 = t22 * t15 + t24 * t18;
t16 = sin(qJ(1));
t19 = cos(qJ(1));
t26 = t19 * t18;
t5 = t11 * t16 - t12 * t26;
t6 = -t11 * t26 - t12 * t16;
t31 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t27 = t16 * t18;
t7 = -t11 * t19 - t12 * t27;
t8 = -t11 * t27 + t12 * t19;
t30 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t17 = cos(qJ(5));
t29 = pkin(5) * t17;
t28 = r_i_i_C(2) * t11;
t25 = pkin(1) + pkin(7) + pkin(4) + t29;
t21 = t23 * t18 + qJ(2) + t33;
t9 = t15 * t12 * r_i_i_C(1);
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t25 * t16 + t21 * t19, t16, t32 * t16, -t27 (-t14 * t19 - t17 * t27) * pkin(5) + t30, t30; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t21 * t16 + t25 * t19, -t19, -t32 * t19, t26 (-t14 * t16 + t17 * t26) * pkin(5) + t31, t31; 0, 0, t22 * t18 - t33, t15, t9 + (-t28 + t29) * t15, -t15 * t28 + t9;];
Ja_transl  = t1;
