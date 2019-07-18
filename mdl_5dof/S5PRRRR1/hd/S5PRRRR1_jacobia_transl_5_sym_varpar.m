% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5PRRRR1
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% Ja_transl [3x5]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobia_transl_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobia_transl_5_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:29:13
% EndTime: 2019-07-18 13:29:14
% DurationCPUTime: 0.12s
% Computational Cost: add. (62->22), mult. (91->36), div. (0->0), fcn. (99->8), ass. (0->27)
t12 = qJ(3) + qJ(4);
t10 = sin(t12);
t11 = cos(t12);
t13 = sin(qJ(5));
t31 = r_i_i_C(2) * t13;
t35 = r_i_i_C(3) * t11 + t10 * t31;
t15 = sin(qJ(2));
t34 = t35 * t15;
t18 = cos(qJ(2));
t33 = t35 * t18;
t16 = cos(qJ(5));
t32 = r_i_i_C(1) * t16;
t29 = t10 * r_i_i_C(3);
t28 = cos(qJ(3)) * pkin(2);
t27 = t15 * t13;
t26 = t15 * t16;
t25 = t18 * t13;
t24 = t18 * t16;
t23 = t10 * t32;
t21 = t28 + t29;
t20 = -sin(qJ(3)) * pkin(2) - t23;
t19 = -t29 + (t31 - t32) * t11;
t4 = t11 * t24 + t27;
t3 = -t11 * t25 + t26;
t2 = -t11 * t26 + t25;
t1 = t11 * t27 + t24;
t5 = [0, t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t21 * t15, t20 * t18 + t33, -t18 * t23 + t33, t3 * r_i_i_C(1) - t4 * r_i_i_C(2); 0, 0, t19 - t28, t19, (r_i_i_C(1) * t13 + r_i_i_C(2) * t16) * t10; 1, t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t21 * t18, t20 * t15 + t34, -t15 * t23 + t34, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2);];
Ja_transl  = t5;
