% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:26:59
% EndTime: 2019-02-26 22:26:59
% DurationCPUTime: 0.12s
% Computational Cost: add. (167->37), mult. (150->51), div. (0->0), fcn. (165->10), ass. (0->30)
t23 = cos(qJ(2));
t21 = sin(qJ(2));
t32 = r_i_i_C(3) + qJ(5) + pkin(9) + pkin(8);
t28 = t32 * t21;
t20 = qJ(3) + qJ(4);
t17 = cos(t20);
t11 = pkin(4) * t17 + cos(qJ(3)) * pkin(3);
t9 = pkin(2) + t11;
t37 = t23 * t9 + pkin(1) + t28;
t15 = pkin(10) + t20;
t12 = sin(t15);
t13 = cos(t15);
t24 = cos(qJ(1));
t22 = sin(qJ(1));
t31 = t22 * t23;
t5 = t12 * t31 + t13 * t24;
t6 = t12 * t24 - t13 * t31;
t36 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t30 = t23 * t24;
t7 = -t12 * t30 + t22 * t13;
t8 = t22 * t12 + t13 * t30;
t35 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t16 = sin(t20);
t34 = pkin(4) * t16;
t10 = t34 + sin(qJ(3)) * pkin(3);
t33 = pkin(7) + t10;
t27 = -r_i_i_C(1) * t12 - r_i_i_C(2) * t13;
t26 = r_i_i_C(1) * t13 - r_i_i_C(2) * t12 + t9;
t25 = -t26 * t21 + t32 * t23;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t37 * t22 + t33 * t24, t25 * t24, -t10 * t30 + t22 * t11 + t35 (-t16 * t30 + t17 * t22) * pkin(4) + t35, t24 * t21, 0; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t33 * t22 + t37 * t24, t25 * t22, -t10 * t31 - t11 * t24 + t36 (-t16 * t31 - t17 * t24) * pkin(4) + t36, t22 * t21, 0; 0, t26 * t23 + t28 (-t10 + t27) * t21 (t27 - t34) * t21, -t23, 0;];
Ja_transl  = t1;
