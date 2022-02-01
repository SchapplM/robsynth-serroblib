% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPRP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(7);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (26->10), mult. (32->14), div. (0->0), fcn. (20->6), ass. (0->6)
	t19 = r_i_i_C(3) + qJ(3);
	t18 = -r_i_i_C(1) * cos(pkin(8)) + r_i_i_C(2) * sin(pkin(8)) - pkin(2);
	t15 = qJ(1) + pkin(7);
	t14 = cos(t15);
	t13 = sin(t15);
	t1 = [t14 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t19 * t13 + t18 * t14) * qJD(1), 0, qJD(1) * t14, 0, 0; t13 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t19 * t14 + t18 * t13) * qJD(1), 0, qJD(1) * t13, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:05
	% EndTime: 2022-01-23 09:13:05
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (89->25), mult. (128->44), div. (0->0), fcn. (106->8), ass. (0->19)
	t144 = cos(pkin(8));
	t145 = sin(qJ(4));
	t153 = t144 * t145;
	t146 = cos(qJ(4));
	t152 = t144 * t146;
	t142 = qJ(1) + pkin(7);
	t140 = sin(t142);
	t141 = cos(t142);
	t151 = -t140 * t145 - t141 * t152;
	t150 = -t140 * t146 + t141 * t153;
	t149 = t140 * t152 - t141 * t145;
	t148 = t140 * t153 + t141 * t146;
	t143 = sin(pkin(8));
	t147 = -pkin(3) * t144 - pkin(2) + (-pkin(6) - r_i_i_C(3)) * t143;
	t139 = t151 * qJD(1) + t148 * qJD(4);
	t138 = t150 * qJD(1) + t149 * qJD(4);
	t137 = t149 * qJD(1) + t150 * qJD(4);
	t136 = t148 * qJD(1) + t151 * qJD(4);
	t1 = [t139 * r_i_i_C(1) + t138 * r_i_i_C(2) + t141 * qJD(3) + (-cos(qJ(1)) * pkin(1) - qJ(3) * t140 + t147 * t141) * qJD(1), 0, qJD(1) * t141, t136 * r_i_i_C(1) + t137 * r_i_i_C(2), 0; -t137 * r_i_i_C(1) + t136 * r_i_i_C(2) + t140 * qJD(3) + (-sin(qJ(1)) * pkin(1) + qJ(3) * t141 + t147 * t140) * qJD(1), 0, qJD(1) * t140, -t138 * r_i_i_C(1) + t139 * r_i_i_C(2), 0; 0, 0, 0, (-r_i_i_C(1) * t146 + r_i_i_C(2) * t145) * t143 * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:05
	% EndTime: 2022-01-23 09:13:05
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (129->34), mult. (189->54), div. (0->0), fcn. (152->8), ass. (0->25)
	t164 = pkin(4) + r_i_i_C(1);
	t146 = qJ(1) + pkin(7);
	t144 = sin(t146);
	t145 = cos(t146);
	t151 = cos(qJ(4));
	t148 = cos(pkin(8));
	t150 = sin(qJ(4));
	t162 = t148 * t150;
	t156 = -t144 * t151 + t145 * t162;
	t163 = t156 * qJD(4);
	t161 = t148 * t151;
	t155 = t144 * t161 - t145 * t150;
	t141 = t156 * qJD(1) + t155 * qJD(4);
	t147 = sin(pkin(8));
	t160 = qJD(1) * t147;
	t159 = qJD(5) * t147;
	t158 = pkin(4) * t150 + qJ(3);
	t157 = -t144 * t150 - t145 * t161;
	t154 = t144 * t162 + t145 * t151;
	t153 = -(pkin(4) * t151 + pkin(3)) * t148 - pkin(2) + (-r_i_i_C(3) - qJ(5) - pkin(6)) * t147;
	t152 = t154 * qJD(4);
	t139 = t154 * qJD(1) + t157 * qJD(4);
	t142 = t157 * qJD(1) + t152;
	t140 = t155 * qJD(1) + t163;
	t1 = [-t144 * t159 + t142 * r_i_i_C(1) + t141 * r_i_i_C(2) + t145 * qJD(3) + pkin(4) * t152 + (-cos(qJ(1)) * pkin(1) - t158 * t144 + t153 * t145) * qJD(1), 0, qJD(1) * t145, t140 * r_i_i_C(2) + t164 * t139, -t144 * t160; t145 * t159 - t140 * r_i_i_C(1) + t139 * r_i_i_C(2) + t144 * qJD(3) - pkin(4) * t163 + (-sin(qJ(1)) * pkin(1) + t158 * t145 + t153 * t144) * qJD(1), 0, qJD(1) * t144, t142 * r_i_i_C(2) - t164 * t141, t145 * t160; 0, 0, 0, (r_i_i_C(2) * t150 - t164 * t151) * t147 * qJD(4), 0;];
	JaD_transl = t1;
end