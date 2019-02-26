% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRP2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:42
% EndTime: 2019-02-26 20:15:42
% DurationCPUTime: 0.18s
% Computational Cost: add. (240->41), mult. (389->71), div. (0->0), fcn. (485->12), ass. (0->36)
t56 = pkin(10) + r_i_i_C(3);
t29 = sin(qJ(5));
t32 = cos(qJ(5));
t55 = t32 * r_i_i_C(1) - t29 * r_i_i_C(2) + pkin(4);
t25 = qJ(3) + qJ(4);
t23 = sin(t25);
t24 = cos(t25);
t33 = cos(qJ(3));
t54 = t33 * pkin(3) + t56 * t23 + t55 * t24 + pkin(2);
t26 = sin(pkin(11));
t27 = sin(pkin(6));
t50 = t26 * t27;
t31 = sin(qJ(2));
t49 = t26 * t31;
t34 = cos(qJ(2));
t48 = t26 * t34;
t47 = t27 * t31;
t46 = t27 * t34;
t45 = cos(pkin(11));
t44 = t27 * t45;
t43 = t45 * t31;
t42 = t45 * t34;
t40 = t29 * r_i_i_C(1) + t32 * r_i_i_C(2) + pkin(8) + pkin(9);
t28 = cos(pkin(6));
t17 = t28 * t43 + t48;
t8 = t17 * t24 - t23 * t44;
t39 = t56 * t8 + t55 * (-t17 * t23 - t24 * t44);
t19 = -t28 * t49 + t42;
t10 = t19 * t24 + t23 * t50;
t38 = t56 * t10 + t55 * (-t19 * t23 + t24 * t50);
t15 = t28 * t23 + t24 * t47;
t37 = t56 * t15 + t55 * (-t23 * t47 + t28 * t24);
t30 = sin(qJ(3));
t18 = t28 * t48 + t43;
t16 = -t28 * t42 + t49;
t1 = [0, -t18 * t54 + t40 * t19 (-t19 * t30 + t33 * t50) * pkin(3) + t38, t38 (-t10 * t29 + t18 * t32) * r_i_i_C(1) + (-t10 * t32 - t18 * t29) * r_i_i_C(2), 0; 0, -t16 * t54 + t40 * t17 (-t17 * t30 - t33 * t44) * pkin(3) + t39, t39 (t16 * t32 - t8 * t29) * r_i_i_C(1) + (-t16 * t29 - t8 * t32) * r_i_i_C(2), 0; 1 (t40 * t31 + t34 * t54) * t27 (t28 * t33 - t30 * t47) * pkin(3) + t37, t37 (-t15 * t29 - t32 * t46) * r_i_i_C(1) + (-t15 * t32 + t29 * t46) * r_i_i_C(2), 0;];
Ja_transl  = t1;
