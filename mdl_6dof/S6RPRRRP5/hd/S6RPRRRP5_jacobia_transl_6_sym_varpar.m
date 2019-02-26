% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP5
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
% Datum: 2019-02-26 21:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:10:18
% EndTime: 2019-02-26 21:10:19
% DurationCPUTime: 0.12s
% Computational Cost: add. (238->33), mult. (198->40), div. (0->0), fcn. (221->9), ass. (0->31)
t26 = sin(qJ(5));
t28 = cos(qJ(5));
t35 = r_i_i_C(3) + qJ(6);
t46 = pkin(5) + r_i_i_C(1);
t52 = t35 * t26 + t46 * t28;
t51 = pkin(9) + r_i_i_C(2);
t25 = pkin(10) + qJ(3);
t23 = qJ(4) + t25;
t19 = cos(t23);
t49 = t19 * t51;
t18 = sin(t23);
t48 = t19 * pkin(4) + t51 * t18;
t17 = pkin(3) * cos(t25);
t47 = t17 + cos(pkin(10)) * pkin(2) + pkin(1) + t48;
t27 = sin(qJ(1));
t44 = t27 * t49;
t39 = t27 * t26;
t38 = t27 * t28;
t29 = cos(qJ(1));
t37 = t29 * t26;
t36 = t29 * t28;
t34 = t29 * t49;
t32 = t52 * t19 + t48;
t31 = (-pkin(4) - t52) * t18;
t30 = -pkin(3) * sin(t25) + t31;
t24 = -pkin(8) - pkin(7) - qJ(2);
t4 = t19 * t36 + t39;
t3 = t19 * t37 - t38;
t2 = t19 * t38 - t37;
t1 = t19 * t39 + t36;
t5 = [-t35 * t1 - t46 * t2 - t29 * t24 - t47 * t27, t27, t30 * t29 + t34, t29 * t31 + t34, -t46 * t3 + t35 * t4, t3; -t27 * t24 + t47 * t29 + t35 * t3 + t46 * t4, -t29, t30 * t27 + t44, t27 * t31 + t44, -t46 * t1 + t35 * t2, t1; 0, 0, t17 + t32, t32 (-t46 * t26 + t35 * t28) * t18, t18 * t26;];
Ja_transl  = t5;
