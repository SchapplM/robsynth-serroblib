% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR14
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR14_jacobia_transl_3_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_transl_3_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobia_transl_3_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_transl_3_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:20
% EndTime: 2018-12-10 18:38:20
% DurationCPUTime: 0.18s
% Computational Cost: add. (231->54), mult. (251->75), div. (0->0), fcn. (251->18), ass. (0->39)
t45 = r_i_i_C(3) + qJ(3);
t26 = sin(pkin(6));
t30 = sin(qJ(1));
t44 = t26 * t30;
t32 = cos(qJ(1));
t43 = t26 * t32;
t42 = pkin(6) - qJ(2);
t41 = pkin(6) + qJ(2);
t25 = sin(pkin(7));
t40 = t45 * t25;
t39 = cos(t42);
t38 = sin(t42);
t37 = cos(t41) / 0.2e1;
t36 = sin(t41) / 0.2e1;
t14 = t36 - t38 / 0.2e1;
t31 = cos(qJ(2));
t35 = -t32 * t14 - t30 * t31;
t7 = t30 * t14 - t32 * t31;
t29 = sin(qJ(2));
t33 = t39 / 0.2e1 + t37;
t2 = t30 * t29 - t32 * t33;
t28 = cos(pkin(7));
t34 = -t2 * t25 + t28 * t43;
t5 = -t32 * t29 - t30 * t33;
t1 = -t5 * t25 + t28 * t44;
t27 = cos(pkin(14));
t24 = sin(pkin(14));
t23 = pkin(7) - pkin(14);
t22 = pkin(7) + pkin(14);
t21 = cos(t22);
t20 = sin(t23);
t17 = cos(t23) / 0.2e1;
t16 = sin(t22) / 0.2e1;
t13 = t36 + t38 / 0.2e1;
t12 = t17 - t21 / 0.2e1;
t11 = t17 + t21 / 0.2e1;
t10 = t16 - t20 / 0.2e1;
t9 = t16 + t20 / 0.2e1;
t3 = [(t2 * t10 + t12 * t43 + t27 * t35) * r_i_i_C(1) + (t2 * t11 - t24 * t35 + t9 * t43) * r_i_i_C(2) + t35 * pkin(2) - t30 * pkin(1) + pkin(10) * t43 + t45 * t34 (t7 * t10 + t5 * t27) * r_i_i_C(1) + (t7 * t11 - t5 * t24) * r_i_i_C(2) + t5 * pkin(2) - t7 * t40, t1, 0, 0, 0; (t5 * t10 + t12 * t44 - t27 * t7) * r_i_i_C(1) + (t5 * t11 + t24 * t7 + t9 * t44) * r_i_i_C(2) - t7 * pkin(2) + t32 * pkin(1) + pkin(10) * t44 + t45 * t1 (t10 * t35 - t2 * t27) * r_i_i_C(1) + (t11 * t35 + t2 * t24) * r_i_i_C(2) - t2 * pkin(2) - t35 * t40, -t34, 0, 0, 0; 0 (t27 * r_i_i_C(1) - t24 * r_i_i_C(2) + pkin(2)) * t13 + (t10 * r_i_i_C(1) + t11 * r_i_i_C(2) - t40) * (t37 - t39 / 0.2e1) -t13 * t25 + cos(pkin(6)) * t28, 0, 0, 0;];
Ja_transl  = t3;
