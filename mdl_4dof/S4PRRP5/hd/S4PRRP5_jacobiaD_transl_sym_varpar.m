% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4PRRP5
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4PRRP5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRRP5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRP5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_jacobiaD_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:29:29
	% EndTime: 2019-12-31 16:29:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:29:29
	% EndTime: 2019-12-31 16:29:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:29:29
	% EndTime: 2019-12-31 16:29:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(6)) * t14, 0, 0; 0, sin(pkin(6)) * t14, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:29:30
	% EndTime: 2019-12-31 16:29:30
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (31->19), mult. (116->43), div. (0->0), fcn. (87->6), ass. (0->17)
	t151 = pkin(5) + r_i_i_C(3);
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
	t137 = cos(pkin(6));
	t136 = sin(pkin(6));
	t1 = [0, t142 * t137, t137 * t143 + ((-t136 * t138 - t137 * t149) * r_i_i_C(1) + (-t136 * t140 + t137 * t150) * r_i_i_C(2)) * qJD(3), 0; 0, t142 * t136, t136 * t143 + ((-t136 * t149 + t137 * t138) * r_i_i_C(1) + (t136 * t150 + t137 * t140) * r_i_i_C(2)) * qJD(3), 0; 0, -t146 * t141 * qJD(3) + (t145 * t139 + t151 * t141) * qJD(2), (t138 * t147 - t140 * t148) * r_i_i_C(2) + (-t138 * t148 - t140 * t147) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:29:30
	% EndTime: 2019-12-31 16:29:30
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (51->20), mult. (168->44), div. (0->0), fcn. (128->6), ass. (0->19)
	t145 = sin(qJ(3));
	t147 = cos(qJ(3));
	t162 = pkin(3) + r_i_i_C(1);
	t164 = r_i_i_C(2) * t145 - t162 * t147;
	t163 = r_i_i_C(2) * t147 + t162 * t145;
	t159 = r_i_i_C(3) + qJ(4) + pkin(5);
	t148 = cos(qJ(2));
	t158 = t145 * t148;
	t157 = t147 * t148;
	t156 = qJD(2) * t148;
	t155 = -pkin(2) + t164;
	t146 = sin(qJ(2));
	t152 = t146 * t163;
	t151 = t163 * t148;
	t150 = qJD(2) * t152;
	t149 = qJD(4) * t148 + qJD(3) * t152 + (-t159 * t146 + t155 * t148) * qJD(2);
	t143 = cos(pkin(6));
	t142 = sin(pkin(6));
	t1 = [0, t149 * t143, t143 * t150 + ((-t142 * t147 + t143 * t158) * r_i_i_C(2) + t162 * (-t142 * t145 - t143 * t157)) * qJD(3), t143 * t156; 0, t149 * t142, t142 * t150 + ((t142 * t158 + t143 * t147) * r_i_i_C(2) + t162 * (-t142 * t157 + t143 * t145)) * qJD(3), t142 * t156; 0, t146 * qJD(4) - qJD(3) * t151 + (t155 * t146 + t159 * t148) * qJD(2), t164 * t146 * qJD(3) - qJD(2) * t151, qJD(2) * t146;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end