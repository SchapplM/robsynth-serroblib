% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:51:42
% EndTime: 2019-02-26 20:51:42
% DurationCPUTime: 0.14s
% Computational Cost: add. (200->36), mult. (135->47), div. (0->0), fcn. (152->11), ass. (0->32)
t21 = pkin(10) + qJ(3);
t17 = cos(t21);
t15 = sin(t21);
t36 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(4);
t29 = t36 * t15;
t20 = pkin(11) + qJ(5);
t16 = cos(t20);
t9 = pkin(5) * t16 + cos(pkin(11)) * pkin(4) + pkin(3);
t40 = t29 + t17 * t9 + cos(pkin(10)) * pkin(2) + pkin(1);
t18 = qJ(6) + t20;
t12 = cos(t18);
t24 = cos(qJ(1));
t31 = t24 * t12;
t11 = sin(t18);
t23 = sin(qJ(1));
t34 = t23 * t11;
t5 = t17 * t34 + t31;
t32 = t24 * t11;
t33 = t23 * t12;
t6 = -t17 * t33 + t32;
t39 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t17 * t32 + t33;
t8 = t17 * t31 + t34;
t38 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t14 = sin(t20);
t37 = pkin(5) * t14;
t35 = t14 * t17;
t30 = t37 + sin(pkin(11)) * pkin(4) + pkin(7) + qJ(2);
t27 = -r_i_i_C(1) * t11 - r_i_i_C(2) * t12;
t26 = r_i_i_C(1) * t12 - r_i_i_C(2) * t11 + t9;
t25 = -t26 * t15 + t36 * t17;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t40 * t23 + t30 * t24, t23, t25 * t24, t24 * t15 (t16 * t23 - t24 * t35) * pkin(5) + t38, t38; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t30 * t23 + t40 * t24, -t24, t25 * t23, t23 * t15 (-t16 * t24 - t23 * t35) * pkin(5) + t39, t39; 0, 0, t26 * t17 + t29, -t17 (t27 - t37) * t15, t27 * t15;];
Ja_transl  = t1;
