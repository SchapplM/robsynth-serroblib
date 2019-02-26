% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:03
% EndTime: 2019-02-26 19:47:03
% DurationCPUTime: 0.13s
% Computational Cost: add. (165->31), mult. (431->56), div. (0->0), fcn. (566->12), ass. (0->30)
t25 = sin(pkin(10));
t26 = sin(pkin(6));
t43 = t25 * t26;
t29 = cos(pkin(10));
t42 = t29 * t26;
t30 = cos(pkin(6));
t34 = cos(qJ(2));
t41 = t30 * t34;
t40 = r_i_i_C(3) + qJ(5);
t24 = sin(pkin(11));
t28 = cos(pkin(11));
t32 = sin(qJ(2));
t38 = t34 * t24 + t32 * t28;
t17 = t38 * t30;
t18 = t32 * t24 - t34 * t28;
t7 = t29 * t17 - t25 * t18;
t39 = t25 * t17 + t29 * t18;
t23 = sin(pkin(12));
t27 = cos(pkin(12));
t37 = r_i_i_C(1) * t27 - r_i_i_C(2) * t23 + pkin(4);
t36 = t23 * r_i_i_C(1) + t27 * r_i_i_C(2) + pkin(8);
t31 = sin(qJ(4));
t33 = cos(qJ(4));
t35 = t40 * t31 + t37 * t33 + pkin(3);
t16 = t18 * t30;
t15 = t38 * t26;
t11 = t15 * t31 - t30 * t33;
t3 = -t31 * t39 - t33 * t43;
t1 = t7 * t31 + t33 * t42;
t2 = [0 (-t25 * t41 - t29 * t32) * pkin(2) - t36 * t39 + t35 * (t25 * t16 - t29 * t38) t43, t40 * (t31 * t43 - t33 * t39) - t37 * t3, t3, 0; 0 (-t25 * t32 + t29 * t41) * pkin(2) + t36 * t7 + t35 * (-t29 * t16 - t25 * t38) -t42, t40 * (-t31 * t42 + t7 * t33) - t37 * t1, t1, 0; 1, t36 * t15 + (pkin(2) * t34 - t18 * t35) * t26, t30, t40 * (t15 * t33 + t30 * t31) - t37 * t11, t11, 0;];
Ja_transl  = t2;
