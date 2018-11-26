% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRRRR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
%
% Output:
% Ja_transl [3x5]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRR1_jacobia_transl_5_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_jacobia_transl_5_floatb_twist_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR1_jacobia_transl_5_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_jacobia_transl_5_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 14:52:29
% EndTime: 2018-11-16 14:52:29
% DurationCPUTime: 0.14s
% Computational Cost: add. (200->32), mult. (160->41), div. (0->0), fcn. (168->10), ass. (0->32)
t19 = qJ(2) + qJ(3);
t18 = qJ(4) + t19;
t14 = sin(t18);
t15 = cos(t18);
t20 = sin(qJ(5));
t42 = r_i_i_C(2) * t20;
t49 = pkin(6) + r_i_i_C(3);
t50 = t14 * t42 + t15 * t49;
t32 = t49 * t14;
t43 = pkin(3) * cos(t19);
t45 = cos(qJ(2)) * pkin(2);
t48 = pkin(4) * t15 + pkin(1) + t32 + t43 + t45;
t22 = cos(qJ(5));
t31 = -r_i_i_C(1) * t22 - pkin(4);
t29 = t31 * t14;
t27 = t29 - pkin(3) * sin(t19);
t24 = cos(qJ(1));
t39 = t20 * t24;
t21 = sin(qJ(1));
t38 = t21 * t20;
t37 = t21 * t22;
t36 = t22 * t24;
t35 = t50 * t21;
t33 = t50 * t24;
t28 = -sin(qJ(2)) * pkin(2) + t27;
t26 = -t32 + (t31 + t42) * t15;
t25 = t26 - t43;
t4 = t15 * t36 - t38;
t3 = -t15 * t39 - t37;
t2 = -t15 * t37 - t39;
t1 = t15 * t38 - t36;
t5 = [r_i_i_C(1) * t2 + r_i_i_C(2) * t1 - t48 * t21, t28 * t24 + t33, t27 * t24 + t33, t24 * t29 + t33, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t48 * t24, t28 * t21 + t35, t27 * t21 + t35, t21 * t29 + t35, -t1 * r_i_i_C(1) + r_i_i_C(2) * t2; 0, t25 - t45, t25, t26 (r_i_i_C(1) * t20 + r_i_i_C(2) * t22) * t14;];
Ja_transl  = t5;
