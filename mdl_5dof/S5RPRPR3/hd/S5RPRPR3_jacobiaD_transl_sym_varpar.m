% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR3
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:42
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (r_i_i_C(1) * t5 + r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t6 + r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t9 = qJ(1) + pkin(8);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; (sin(qJ(1)) * pkin(1) + r_i_i_C(1) * t7 + r_i_i_C(2) * t8) * qJD(1), 0, 0, 0, 0; (-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t8 + r_i_i_C(2) * t7) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (34->9), mult. (24->12), div. (0->0), fcn. (12->6), ass. (0->8)
	t18 = qJ(1) + pkin(8);
	t16 = qJ(3) + t18;
	t17 = qJD(1) + qJD(3);
	t22 = sin(t16) * t17;
	t21 = cos(t16) * t17;
	t20 = r_i_i_C(1) * t22 + r_i_i_C(2) * t21;
	t19 = -r_i_i_C(1) * t21 + r_i_i_C(2) * t22;
	t1 = [0, 0, 0, 0, 0; (sin(t18) * pkin(2) + sin(qJ(1)) * pkin(1)) * qJD(1) + t20, 0, t20, 0, 0; (-cos(t18) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t19, 0, t19, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (101->16), mult. (62->20), div. (0->0), fcn. (38->8), ass. (0->13)
	t70 = pkin(3) + r_i_i_C(1) * cos(pkin(9));
	t68 = r_i_i_C(2) * sin(pkin(9));
	t60 = qJ(1) + pkin(8);
	t58 = qJ(3) + t60;
	t56 = sin(t58);
	t59 = qJD(1) + qJD(3);
	t67 = t59 * t56;
	t57 = cos(t58);
	t66 = t59 * t57;
	t65 = -r_i_i_C(3) - qJ(4);
	t64 = t66 * t68 + t57 * qJD(4) + (t65 * t56 - t57 * t70) * t59;
	t63 = -t56 * qJD(4) + (-t56 * t68 + t65 * t57) * t59 + t70 * t67;
	t1 = [0, 0, 0, 0, 0; (sin(t60) * pkin(2) + sin(qJ(1)) * pkin(1)) * qJD(1) + t63, 0, t63, -t67, 0; (-cos(t60) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t64, 0, t64, t66, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (268->30), mult. (204->47), div. (0->0), fcn. (166->10), ass. (0->25)
	t145 = sin(pkin(9));
	t146 = cos(pkin(9));
	t163 = pkin(4) * t146 + (pkin(7) + r_i_i_C(3)) * t145 + pkin(3);
	t144 = qJ(1) + pkin(8);
	t142 = qJ(3) + t144;
	t140 = sin(t142);
	t143 = qJD(1) + qJD(3);
	t159 = t143 * t140;
	t141 = cos(t142);
	t158 = t143 * t141;
	t147 = sin(qJ(5));
	t157 = t146 * t147;
	t148 = cos(qJ(5));
	t156 = t146 * t148;
	t154 = t140 * t147 + t141 * t156;
	t153 = -t140 * t148 + t141 * t157;
	t152 = t140 * t156 - t141 * t147;
	t151 = t140 * t157 + t141 * t148;
	t131 = -t154 * qJD(5) + t151 * t143;
	t132 = t153 * qJD(5) + t152 * t143;
	t150 = t132 * r_i_i_C(1) - t131 * r_i_i_C(2) - qJ(4) * t158 - t140 * qJD(4) + t163 * t159;
	t133 = t152 * qJD(5) + t153 * t143;
	t134 = -t151 * qJD(5) + t154 * t143;
	t149 = t133 * r_i_i_C(2) - t134 * r_i_i_C(1) + t141 * qJD(4) + (-qJ(4) * t140 - t163 * t141) * t143;
	t1 = [0, 0, 0, 0, (-r_i_i_C(1) * t148 + r_i_i_C(2) * t147) * t145 * qJD(5); (sin(t144) * pkin(2) + sin(qJ(1)) * pkin(1)) * qJD(1) + t150, 0, t150, -t159, t133 * r_i_i_C(1) + t134 * r_i_i_C(2); (-cos(t144) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t149, 0, t149, t158, t131 * r_i_i_C(1) + t132 * r_i_i_C(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end