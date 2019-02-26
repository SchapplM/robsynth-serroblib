% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR4_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR4_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:10
% EndTime: 2019-02-26 21:30:10
% DurationCPUTime: 0.10s
% Computational Cost: add. (92->29), mult. (216->43), div. (0->0), fcn. (281->8), ass. (0->24)
t32 = pkin(3) - r_i_i_C(2);
t22 = cos(qJ(2));
t31 = pkin(2) * t22;
t19 = cos(pkin(6));
t30 = t19 * t22;
t29 = r_i_i_C(3) + qJ(4);
t17 = sin(pkin(6));
t20 = sin(qJ(2));
t28 = -pkin(2) * t19 * t20 + (r_i_i_C(1) + pkin(8) + qJ(3)) * t17;
t16 = sin(pkin(11));
t18 = cos(pkin(11));
t25 = t16 * t22 + t18 * t20;
t10 = t25 * t19;
t12 = t16 * t20 - t22 * t18;
t21 = sin(qJ(1));
t23 = cos(qJ(1));
t27 = t10 * t23 - t21 * t12;
t26 = t21 * t10 + t23 * t12;
t24 = t12 * t19;
t15 = pkin(1) + t31;
t8 = t12 * t17;
t5 = t21 * t24 - t23 * t25;
t2 = -t21 * t25 - t23 * t24;
t1 = [-t21 * t15 + t29 * t2 + t28 * t23 - t32 * t27, t32 * t5 - t29 * t26 + (-t20 * t23 - t21 * t30) * pkin(2), t21 * t17, -t5, 0, 0; t15 * t23 + t28 * t21 - t32 * t26 - t29 * t5, t32 * t2 + t29 * t27 + (-t20 * t21 + t23 * t30) * pkin(2), -t23 * t17, -t2, 0, 0; 0, -t32 * t8 + (t25 * t29 + t31) * t17, t19, t8, 0, 0;];
Ja_transl  = t1;
