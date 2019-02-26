% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:51
% EndTime: 2019-02-26 21:07:51
% DurationCPUTime: 0.15s
% Computational Cost: add. (235->33), mult. (198->44), div. (0->0), fcn. (219->10), ass. (0->28)
t42 = pkin(5) + r_i_i_C(1);
t46 = pkin(9) + r_i_i_C(2);
t34 = r_i_i_C(3) + qJ(6);
t25 = qJ(3) + qJ(4);
t22 = cos(t25);
t45 = t22 * t46;
t21 = sin(t25);
t44 = t22 * pkin(4) + t46 * t21;
t23 = cos(qJ(3)) * pkin(3);
t43 = pkin(2) + t23 + t44;
t24 = qJ(1) + pkin(10);
t19 = sin(t24);
t40 = t19 * t45;
t20 = cos(t24);
t39 = t20 * t45;
t26 = sin(qJ(5));
t36 = t22 * t26;
t28 = cos(qJ(5));
t35 = t22 * t28;
t32 = t34 * t36 + t42 * t35 + t44;
t31 = (-t34 * t26 - t42 * t28 - pkin(4)) * t21;
t30 = -sin(qJ(3)) * pkin(3) + t31;
t29 = -pkin(8) - pkin(7);
t4 = t19 * t26 + t20 * t35;
t3 = -t19 * t28 + t20 * t36;
t2 = t19 * t35 - t20 * t26;
t1 = t19 * t36 + t20 * t28;
t5 = [-sin(qJ(1)) * pkin(1) - t20 * t29 - t42 * t2 - t34 * t1 - t43 * t19, 0, t20 * t30 + t39, t20 * t31 + t39, -t42 * t3 + t34 * t4, t3; cos(qJ(1)) * pkin(1) - t19 * t29 + t42 * t4 + t34 * t3 + t43 * t20, 0, t19 * t30 + t40, t19 * t31 + t40, -t42 * t1 + t34 * t2, t1; 0, 1, t23 + t32, t32 (-t42 * t26 + t34 * t28) * t21, t21 * t26;];
Ja_transl  = t5;
