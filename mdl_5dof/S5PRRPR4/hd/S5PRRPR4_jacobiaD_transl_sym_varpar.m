% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:08
	% EndTime: 2019-12-05 16:24:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:08
	% EndTime: 2019-12-05 16:24:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:08
	% EndTime: 2019-12-05 16:24:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(8)) * t14, 0, 0, 0; 0, sin(pkin(8)) * t14, 0, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:09
	% EndTime: 2019-12-05 16:24:09
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (31->19), mult. (116->43), div. (0->0), fcn. (87->6), ass. (0->17)
	t151 = pkin(6) + r_i_i_C(3);
	t138 = sin(qJ(3));
	t141 = cos(qJ(2));
	t150 = t138 * t141;
	t140 = cos(qJ(3));
	t149 = t140 * t141;
	t148 = qJD(2) * t141;
	t139 = sin(qJ(2));
	t147 = qJD(3) * t139;
	t146 = r_i_i_C(1) * t138 + r_i_i_C(2) * t140;
	t145 = -r_i_i_C(1) * t140 + r_i_i_C(2) * t138 - pkin(2);
	t144 = t146 * t139;
	t143 = qJD(2) * t144;
	t142 = qJD(3) * t144 + (-t151 * t139 + t145 * t141) * qJD(2);
	t137 = cos(pkin(8));
	t136 = sin(pkin(8));
	t1 = [0, t142 * t137, t137 * t143 + ((-t136 * t138 - t137 * t149) * r_i_i_C(1) + (-t136 * t140 + t137 * t150) * r_i_i_C(2)) * qJD(3), 0, 0; 0, t142 * t136, t136 * t143 + ((-t136 * t149 + t137 * t138) * r_i_i_C(1) + (t136 * t150 + t137 * t140) * r_i_i_C(2)) * qJD(3), 0, 0; 0, -t146 * t141 * qJD(3) + (t145 * t139 + t151 * t141) * qJD(2), (t138 * t147 - t140 * t148) * r_i_i_C(2) + (-t138 * t148 - t140 * t147) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:09
	% EndTime: 2019-12-05 16:24:09
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (79->26), mult. (168->53), div. (0->0), fcn. (128->8), ass. (0->22)
	t145 = qJ(3) + pkin(9);
	t143 = sin(t145);
	t144 = cos(t145);
	t151 = cos(qJ(3));
	t166 = -t151 * pkin(3) - r_i_i_C(1) * t144 + r_i_i_C(2) * t143;
	t164 = r_i_i_C(3) + qJ(4) + pkin(6);
	t146 = sin(pkin(8));
	t152 = cos(qJ(2));
	t163 = t146 * t152;
	t147 = cos(pkin(8));
	t162 = t147 * t152;
	t161 = t151 * t152;
	t160 = qJD(2) * t152;
	t158 = -pkin(2) + t166;
	t149 = sin(qJ(3));
	t157 = pkin(3) * t149 + r_i_i_C(1) * t143 + r_i_i_C(2) * t144;
	t156 = t157 * t152;
	t150 = sin(qJ(2));
	t155 = t157 * t150;
	t154 = qJD(2) * t155;
	t153 = qJD(4) * t152 + qJD(3) * t155 + (-t164 * t150 + t158 * t152) * qJD(2);
	t1 = [0, t153 * t147, t147 * t154 + ((-t143 * t146 - t144 * t162) * r_i_i_C(1) + (t143 * t162 - t144 * t146) * r_i_i_C(2) + (-t146 * t149 - t147 * t161) * pkin(3)) * qJD(3), t147 * t160, 0; 0, t153 * t146, t146 * t154 + ((t143 * t147 - t144 * t163) * r_i_i_C(1) + (t143 * t163 + t144 * t147) * r_i_i_C(2) + (-t146 * t161 + t147 * t149) * pkin(3)) * qJD(3), t146 * t160, 0; 0, t150 * qJD(4) - qJD(3) * t156 + (t158 * t150 + t164 * t152) * qJD(2), t166 * t150 * qJD(3) - qJD(2) * t156, qJD(2) * t150, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:09
	% EndTime: 2019-12-05 16:24:09
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (209->35), mult. (238->61), div. (0->0), fcn. (182->10), ass. (0->34)
	t182 = qJ(3) + pkin(9);
	t174 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t182);
	t170 = t174 * qJD(3);
	t181 = qJD(3) + qJD(5);
	t179 = qJ(5) + t182;
	t175 = sin(t179);
	t176 = cos(t179);
	t205 = r_i_i_C(2) * t176;
	t194 = r_i_i_C(1) * t175 + t205;
	t208 = t194 * t181 - t170;
	t207 = r_i_i_C(1) * t176;
	t206 = r_i_i_C(2) * t175;
	t204 = r_i_i_C(3) + pkin(7) + qJ(4) + pkin(6);
	t186 = sin(qJ(2));
	t203 = t181 * t186;
	t188 = cos(qJ(2));
	t202 = t181 * t188;
	t183 = sin(pkin(8));
	t184 = cos(pkin(8));
	t199 = qJD(2) * t186;
	t190 = t181 * t184 + t183 * t199;
	t197 = t183 * t202;
	t201 = (t190 * t175 - t176 * t197) * r_i_i_C(1) + (t175 * t197 + t190 * t176) * r_i_i_C(2);
	t191 = -t181 * t183 + t184 * t199;
	t196 = t184 * t202;
	t200 = (t191 * t175 - t176 * t196) * r_i_i_C(1) + (t175 * t196 + t191 * t176) * r_i_i_C(2);
	t198 = qJD(2) * t188;
	t195 = -cos(qJ(3)) * pkin(3) - pkin(4) * cos(t182);
	t193 = -pkin(2) + t195 + t206 - t207;
	t171 = t195 * qJD(3);
	t192 = t171 * t188 - t174 * t199;
	t189 = qJD(4) * t188 + t208 * t186 + (-t204 * t186 + t193 * t188) * qJD(2);
	t172 = t203 * t206;
	t1 = [0, t189 * t184, t183 * t170 + t192 * t184 + t200, t184 * t198, t200; 0, t189 * t183, -t184 * t170 + t192 * t183 + t201, t183 * t198, t201; 0, t186 * qJD(4) - t208 * t188 + (t193 * t186 + t204 * t188) * qJD(2), t172 + (-t181 * t207 + t171) * t186 + (t174 - t194) * t198, t199, -t198 * t205 + t172 + (-t175 * t198 - t176 * t203) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end