% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR8
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
% Datum: 2019-02-26 20:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:52:46
% EndTime: 2019-02-26 20:52:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (143->33), mult. (135->44), div. (0->0), fcn. (151->10), ass. (0->31)
t14 = qJ(3) + pkin(10);
t11 = cos(t14);
t33 = r_i_i_C(3) + pkin(9) + pkin(8);
t40 = t33 * t11 - sin(qJ(3)) * pkin(3);
t10 = sin(t14);
t15 = qJ(5) + qJ(6);
t12 = sin(t15);
t13 = cos(t15);
t20 = cos(qJ(5));
t9 = pkin(5) * t20 + pkin(4);
t25 = r_i_i_C(1) * t13 - r_i_i_C(2) * t12 + t9;
t39 = t33 * t10 + t25 * t11 + cos(qJ(3)) * pkin(3);
t22 = cos(qJ(1));
t28 = t13 * t22;
t19 = sin(qJ(1));
t31 = t12 * t19;
t5 = -t10 * t31 + t28;
t29 = t13 * t19;
t30 = t12 * t22;
t6 = t10 * t29 + t30;
t38 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t7 = t10 * t30 + t29;
t8 = t10 * t28 - t31;
t37 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
t17 = sin(qJ(5));
t35 = pkin(5) * t17;
t32 = t10 * t17;
t27 = pkin(1) + qJ(4) + pkin(7) + t35;
t26 = -r_i_i_C(1) * t12 - r_i_i_C(2) * t13;
t24 = t10 * t9 + qJ(2) - t40;
t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) - t27 * t19 + t24 * t22, t19, t39 * t19, t22 (-t19 * t32 + t20 * t22) * pkin(5) + t38, t38; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t24 * t19 + t27 * t22, -t22, -t39 * t22, t19 (t19 * t20 + t22 * t32) * pkin(5) + t37, t37; 0, 0, -t25 * t10 + t40, 0 (t26 - t35) * t11, t26 * t11;];
Ja_transl  = t1;
