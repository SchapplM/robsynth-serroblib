% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:15
% EndTime: 2019-02-26 22:31:16
% DurationCPUTime: 0.12s
% Computational Cost: add. (194->25), mult. (144->27), div. (0->0), fcn. (151->10), ass. (0->23)
t46 = r_i_i_C(3) + qJ(5);
t20 = qJ(2) + qJ(3);
t17 = qJ(4) + t20;
t14 = sin(t17);
t15 = cos(t17);
t22 = cos(pkin(11));
t32 = -r_i_i_C(1) * t22 - pkin(4);
t21 = sin(pkin(11));
t39 = r_i_i_C(2) * t21;
t27 = t46 * t14 + (-t32 - t39) * t15;
t25 = pkin(3) * cos(t20) + t27;
t45 = cos(qJ(2)) * pkin(2) + t25;
t44 = t14 * t39 + t46 * t15;
t31 = t32 * t14;
t26 = t31 - pkin(3) * sin(t20);
t42 = pkin(1) + t45;
t23 = sin(qJ(1));
t35 = t44 * t23;
t24 = cos(qJ(1));
t34 = t44 * t24;
t29 = t21 * r_i_i_C(1) + t22 * r_i_i_C(2) + pkin(7) + pkin(8) + pkin(9);
t28 = -sin(qJ(2)) * pkin(2) + t26;
t1 = [-t42 * t23 + t29 * t24, t28 * t24 + t34, t26 * t24 + t34, t24 * t31 + t34, t24 * t14, 0; t29 * t23 + t42 * t24, t28 * t23 + t35, t26 * t23 + t35, t23 * t31 + t35, t23 * t14, 0; 0, t45, t25, t27, -t15, 0;];
Ja_transl  = t1;
