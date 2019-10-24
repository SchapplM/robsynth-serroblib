% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:29
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRPP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:09
	% EndTime: 2019-10-24 10:29:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:09
	% EndTime: 2019-10-24 10:29:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:09
	% EndTime: 2019-10-24 10:29:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(7)) * t14, 0, 0, 0; 0, sin(pkin(7)) * t14, 0, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:09
	% EndTime: 2019-10-24 10:29:09
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
	t137 = cos(pkin(7));
	t136 = sin(pkin(7));
	t1 = [0, t142 * t137, t137 * t143 + ((-t136 * t138 - t137 * t149) * r_i_i_C(1) + (-t136 * t140 + t137 * t150) * r_i_i_C(2)) * qJD(3), 0, 0; 0, t142 * t136, t136 * t143 + ((-t136 * t149 + t137 * t138) * r_i_i_C(1) + (t136 * t150 + t137 * t140) * r_i_i_C(2)) * qJD(3), 0, 0; 0, -t146 * t141 * qJD(3) + (t145 * t139 + t151 * t141) * qJD(2), (t138 * t147 - t140 * t148) * r_i_i_C(2) + (-t138 * t148 - t140 * t147) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:09
	% EndTime: 2019-10-24 10:29:09
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (79->26), mult. (168->53), div. (0->0), fcn. (128->8), ass. (0->22)
	t145 = qJ(3) + pkin(8);
	t143 = sin(t145);
	t144 = cos(t145);
	t151 = cos(qJ(3));
	t166 = -t151 * pkin(3) - r_i_i_C(1) * t144 + r_i_i_C(2) * t143;
	t164 = r_i_i_C(3) + qJ(4) + pkin(6);
	t146 = sin(pkin(7));
	t152 = cos(qJ(2));
	t163 = t146 * t152;
	t147 = cos(pkin(7));
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
	% StartTime: 2019-10-24 10:29:10
	% EndTime: 2019-10-24 10:29:10
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (164->38), mult. (298->65), div. (0->0), fcn. (240->8), ass. (0->29)
	t193 = qJ(3) + pkin(8);
	t191 = sin(t193);
	t192 = cos(t193);
	t197 = sin(qJ(3));
	t216 = r_i_i_C(3) + qJ(5);
	t220 = pkin(4) + r_i_i_C(1);
	t221 = pkin(3) * t197 + t220 * t191 - t216 * t192;
	t224 = t221 * qJD(3) - (r_i_i_C(2) + qJ(4) + pkin(6)) * qJD(2) - t191 * qJD(5);
	t199 = cos(qJ(3));
	t223 = -t199 * pkin(3) - t216 * t191 - t220 * t192;
	t194 = sin(pkin(7));
	t200 = cos(qJ(2));
	t215 = t194 * t200;
	t195 = cos(pkin(7));
	t214 = t195 * t200;
	t213 = t199 * t200;
	t198 = sin(qJ(2));
	t212 = qJD(2) * t198;
	t211 = qJD(2) * t200;
	t209 = t194 * t212;
	t208 = t195 * t212;
	t207 = t197 * t212;
	t206 = -t191 * t194 - t192 * t214;
	t205 = t191 * t195 - t192 * t215;
	t202 = qJD(4) + (-pkin(2) + t223) * qJD(2);
	t201 = t198 * t224 + t202 * t200;
	t188 = qJD(3) * t206 + t191 * t208;
	t186 = qJD(3) * t205 + t191 * t209;
	t1 = [0, t201 * t195, -t206 * qJD(5) - t216 * (t192 * t208 + (t191 * t214 - t192 * t194) * qJD(3)) + t220 * t188 + (t195 * t207 + (-t194 * t197 - t195 * t213) * qJD(3)) * pkin(3), t195 * t211, -t188; 0, t201 * t194, -t205 * qJD(5) - t216 * (t192 * t209 + (t191 * t215 + t192 * t195) * qJD(3)) + t220 * t186 + (t194 * t207 + (-t194 * t213 + t195 * t197) * qJD(3)) * pkin(3), t194 * t211, -t186; 0, t202 * t198 - t200 * t224, -t221 * t211 + (qJD(3) * t223 + qJD(5) * t192) * t198, t212, qJD(3) * t192 * t198 + t191 * t211;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end