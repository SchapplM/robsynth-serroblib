% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP10_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP10_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:14:19
% EndTime: 2019-02-26 22:14:19
% DurationCPUTime: 0.14s
% Computational Cost: add. (132->41), mult. (334->69), div. (0->0), fcn. (423->10), ass. (0->29)
t20 = sin(qJ(3));
t23 = cos(qJ(3));
t17 = sin(pkin(11));
t19 = cos(pkin(11));
t28 = r_i_i_C(1) * t19 - r_i_i_C(2) * t17 + pkin(3);
t33 = r_i_i_C(3) + qJ(4);
t37 = t33 * t20 + t28 * t23 + pkin(2);
t18 = sin(pkin(6));
t22 = sin(qJ(1));
t36 = t18 * t22;
t35 = t18 * t23;
t25 = cos(qJ(1));
t34 = t18 * t25;
t32 = cos(pkin(6));
t21 = sin(qJ(2));
t24 = cos(qJ(2));
t29 = t25 * t32;
t10 = t21 * t29 + t22 * t24;
t31 = t10 * t23 - t20 * t34;
t30 = t22 * t32;
t27 = t17 * r_i_i_C(1) + t19 * r_i_i_C(2) + pkin(9);
t1 = t10 * t20 + t23 * t34;
t12 = -t21 * t30 + t25 * t24;
t11 = t25 * t21 + t24 * t30;
t9 = t22 * t21 - t24 * t29;
t7 = t18 * t21 * t20 - t32 * t23;
t6 = t12 * t23 + t20 * t36;
t5 = t12 * t20 - t22 * t35;
t2 = [(-t9 * t17 - t19 * t31) * r_i_i_C(1) + (t17 * t31 - t9 * t19) * r_i_i_C(2) - t31 * pkin(3) - t10 * pkin(2) - t9 * pkin(9) - t22 * pkin(1) + pkin(8) * t34 - t33 * t1, -t11 * t37 + t27 * t12, -t28 * t5 + t33 * t6, t5, 0, 0; (t11 * t17 + t6 * t19) * r_i_i_C(1) + (t11 * t19 - t6 * t17) * r_i_i_C(2) + t6 * pkin(3) + t12 * pkin(2) + t11 * pkin(9) + t25 * pkin(1) + pkin(8) * t36 + t33 * t5, t27 * t10 - t37 * t9, -t28 * t1 + t33 * t31, t1, 0, 0; 0 (t27 * t21 + t37 * t24) * t18, t33 * (t32 * t20 + t21 * t35) - t28 * t7, t7, 0, 0;];
Ja_transl  = t2;
