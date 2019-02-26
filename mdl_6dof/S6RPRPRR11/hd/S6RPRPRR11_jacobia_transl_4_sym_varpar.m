% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR11
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR11_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR11_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobia_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:30
% EndTime: 2019-02-26 20:54:31
% DurationCPUTime: 0.14s
% Computational Cost: add. (147->43), mult. (401->73), div. (0->0), fcn. (523->12), ass. (0->39)
t48 = cos(qJ(3));
t22 = sin(pkin(7));
t28 = sin(qJ(3));
t47 = t22 * t28;
t26 = cos(pkin(7));
t46 = t26 * t28;
t21 = sin(pkin(12));
t29 = sin(qJ(1));
t45 = t29 * t21;
t23 = sin(pkin(6));
t44 = t29 * t23;
t25 = cos(pkin(12));
t43 = t29 * t25;
t30 = cos(qJ(1));
t42 = t30 * t21;
t41 = t30 * t23;
t40 = t30 * t25;
t39 = r_i_i_C(3) + qJ(4);
t38 = t23 * qJ(2);
t37 = t22 * t48;
t36 = t26 * t48;
t35 = t23 * t37;
t20 = sin(pkin(13));
t24 = cos(pkin(13));
t34 = -t24 * r_i_i_C(1) + t20 * r_i_i_C(2) - pkin(3);
t27 = cos(pkin(6));
t13 = -t27 * t40 + t45;
t9 = -t13 * t22 + t26 * t41;
t33 = t27 * t43 + t42;
t32 = t33 * t26;
t14 = t27 * t42 + t43;
t31 = -t13 * t46 + t14 * t48 - t41 * t47;
t10 = t33 * t22 + t26 * t44;
t1 = t13 * t36 + t14 * t28 + t30 * t35;
t15 = -t27 * t45 + t40;
t7 = -t27 * t37 + (t21 * t28 - t25 * t36) * t23;
t6 = t15 * t48 + (t22 * t44 - t32) * t28;
t5 = t15 * t28 - t29 * t35 + t48 * t32;
t2 = [(t9 * t20 - t24 * t31) * r_i_i_C(1) + (t20 * t31 + t9 * t24) * r_i_i_C(2) - t31 * pkin(3) - t14 * pkin(2) - t29 * pkin(1) + t30 * t38 - t39 * t1 + t9 * pkin(9), t44, t34 * t5 + t39 * t6, t5, 0, 0; (t10 * t20 + t6 * t24) * r_i_i_C(1) + (t10 * t24 - t6 * t20) * r_i_i_C(2) + t6 * pkin(3) + t15 * pkin(2) + t30 * pkin(1) + t29 * t38 + t39 * t5 + t10 * pkin(9), -t41, t34 * t1 + t39 * t31, t1, 0, 0; 0, t27, t39 * (t27 * t47 + (t48 * t21 + t25 * t46) * t23) + t34 * t7, t7, 0, 0;];
Ja_transl  = t2;
