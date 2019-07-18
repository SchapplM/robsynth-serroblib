% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
%
% Output:
% Ja_transl [3x5]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRR3_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobia_transl_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR3_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobia_transl_4_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:19:38
% EndTime: 2019-07-18 17:19:38
% DurationCPUTime: 0.14s
% Computational Cost: add. (94->25), mult. (119->37), div. (0->0), fcn. (127->8), ass. (0->28)
t17 = qJ(2) + qJ(3);
t15 = sin(t17);
t16 = cos(t17);
t18 = sin(qJ(4));
t37 = r_i_i_C(2) * t18;
t44 = pkin(5) + r_i_i_C(3);
t45 = t15 * t37 + t16 * t44;
t42 = t16 * pkin(2) + t44 * t15;
t39 = cos(qJ(2)) * pkin(1);
t41 = t39 + t42;
t21 = cos(qJ(4));
t38 = r_i_i_C(1) * t21;
t23 = cos(qJ(1));
t34 = t18 * t23;
t20 = sin(qJ(1));
t33 = t20 * t18;
t32 = t20 * t21;
t31 = t21 * t23;
t30 = t45 * t20;
t28 = t45 * t23;
t27 = (-pkin(2) - t38) * t15;
t25 = (-t37 + t38) * t16 + t42;
t24 = -sin(qJ(2)) * pkin(1) + t27;
t4 = t16 * t31 + t33;
t3 = -t16 * t34 + t32;
t2 = -t16 * t32 + t34;
t1 = t16 * t33 + t31;
t5 = [r_i_i_C(1) * t2 + r_i_i_C(2) * t1 - t41 * t20, t24 * t23 + t28, t23 * t27 + t28, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t41 * t23, t24 * t20 + t30, t20 * t27 + t30, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, t25 + t39, t25, (-r_i_i_C(1) * t18 - r_i_i_C(2) * t21) * t15, 0;];
Ja_transl  = t5;
