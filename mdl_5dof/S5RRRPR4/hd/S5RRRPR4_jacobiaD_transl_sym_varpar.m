% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:12:10
	% EndTime: 2019-12-31 21:12:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:12:10
	% EndTime: 2019-12-31 21:12:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:12:10
	% EndTime: 2019-12-31 21:12:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t40 = qJ(1) + qJ(2);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(1) + qJD(2);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [-cos(qJ(1)) * t43 + t42, t42, 0, 0, 0; -sin(qJ(1)) * t43 + t41, t41, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:12:10
	% EndTime: 2019-12-31 21:12:10
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (93->19), mult. (104->33), div. (0->0), fcn. (64->6), ass. (0->19)
	t61 = r_i_i_C(3) + pkin(7);
	t42 = qJ(1) + qJ(2);
	t39 = sin(t42);
	t40 = cos(t42);
	t44 = cos(qJ(3));
	t53 = qJD(3) * t44;
	t41 = qJD(1) + qJD(2);
	t43 = sin(qJ(3));
	t57 = t41 * t43;
	t60 = t39 * t57 - t40 * t53;
	t59 = t39 * t53 + t40 * t57;
	t56 = t41 * t44;
	t55 = pkin(1) * qJD(1);
	t54 = qJD(3) * t43;
	t50 = t39 * t54;
	t47 = t39 * t56 + t40 * t54;
	t46 = r_i_i_C(1) * t50 + ((-r_i_i_C(1) * t44 - pkin(2)) * t40 - t61 * t39) * t41 + t59 * r_i_i_C(2);
	t45 = -t47 * r_i_i_C(1) + t60 * r_i_i_C(2) + (-pkin(2) * t39 + t40 * t61) * t41;
	t1 = [-cos(qJ(1)) * t55 + t46, t46, t60 * r_i_i_C(1) + t47 * r_i_i_C(2), 0, 0; -sin(qJ(1)) * t55 + t45, t45, (-t40 * t56 + t50) * r_i_i_C(2) - t59 * r_i_i_C(1), 0, 0; 0, 0, (-r_i_i_C(1) * t43 - r_i_i_C(2) * t44) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:12:11
	% EndTime: 2019-12-31 21:12:11
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (176->24), mult. (202->34), div. (0->0), fcn. (136->6), ass. (0->26)
	t173 = qJ(1) + qJ(2);
	t171 = cos(t173);
	t172 = qJD(1) + qJD(2);
	t190 = t171 * t172;
	t194 = pkin(3) + r_i_i_C(1);
	t195 = r_i_i_C(2) + pkin(7);
	t193 = r_i_i_C(3) + qJ(4);
	t175 = cos(qJ(3));
	t182 = t193 * t175;
	t174 = sin(qJ(3));
	t185 = t194 * t174;
	t180 = -t185 + t182;
	t192 = pkin(1) * qJD(1);
	t170 = sin(t173);
	t191 = t170 * t172;
	t189 = t172 * t174;
	t188 = qJD(3) * t174;
	t187 = qJD(3) * t175;
	t186 = t174 * qJD(4);
	t183 = t171 * t187;
	t181 = -t193 * t174 - t194 * t175;
	t179 = -pkin(2) + t181;
	t178 = t181 * qJD(3) + qJD(4) * t175;
	t177 = t179 * t191 + t195 * t190 + t193 * t183 + (-qJD(3) * t185 + t186) * t171;
	t176 = t179 * t190 + (-qJD(3) * t182 - t195 * t172 + t194 * t188 - t186) * t170;
	t1 = [-cos(qJ(1)) * t192 + t176, t176, t178 * t171 - t180 * t191, -t170 * t189 + t183, 0; -sin(qJ(1)) * t192 + t177, t177, t178 * t170 + t180 * t190, t170 * t187 + t171 * t189, 0; 0, 0, t180 * qJD(3) + t186, t188, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:12:10
	% EndTime: 2019-12-31 21:12:10
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (394->44), mult. (512->65), div. (0->0), fcn. (435->8), ass. (0->41)
	t93 = qJ(1) + qJ(2);
	t91 = cos(t93);
	t92 = qJD(1) + qJD(2);
	t126 = t91 * t92;
	t94 = sin(qJ(5));
	t95 = sin(qJ(3));
	t96 = cos(qJ(5));
	t97 = cos(qJ(3));
	t107 = t94 * t97 - t95 * t96;
	t128 = pkin(3) + pkin(4);
	t133 = pkin(8) + r_i_i_C(3);
	t90 = sin(t93);
	t127 = t90 * t92;
	t106 = t94 * t95 + t96 * t97;
	t130 = qJD(3) - qJD(5);
	t78 = t130 * t106;
	t118 = t95 * qJD(4);
	t119 = qJD(3) * t97;
	t132 = -qJ(4) * t119 - t118;
	t117 = t128 * t95;
	t131 = -qJ(4) * t97 + t117;
	t120 = qJD(3) * t95;
	t129 = t107 * qJD(5) + t96 * t120;
	t125 = t92 * t95;
	t122 = pkin(1) * qJD(1);
	t114 = t91 * t119;
	t103 = t107 * t92;
	t73 = t90 * t103 + t78 * t91;
	t102 = t106 * t92;
	t74 = t90 * t102 - t94 * t114 + t129 * t91;
	t110 = t73 * r_i_i_C(1) + t74 * r_i_i_C(2);
	t75 = t91 * t103 - t90 * t78;
	t76 = t91 * t102 + (t94 * t119 - t129) * t90;
	t109 = -t75 * r_i_i_C(1) - t76 * r_i_i_C(2);
	t108 = -t130 * t107 * r_i_i_C(1) - t78 * r_i_i_C(2);
	t105 = -qJ(4) * t95 - t128 * t97;
	t101 = -pkin(2) + t105;
	t100 = t105 * qJD(3) + qJD(4) * t97;
	t99 = t101 * t127 + pkin(7) * t126 - t74 * r_i_i_C(1) + t73 * r_i_i_C(2) + (-qJD(3) * t117 - t133 * t92 - t132) * t91;
	t98 = t101 * t126 - t76 * r_i_i_C(1) + t75 * r_i_i_C(2) + t133 * t127 + (-pkin(7) * t92 + t128 * t120 + t132) * t90;
	t1 = [-cos(qJ(1)) * t122 + t98, t98, t100 * t91 + t131 * t127 - t110, -t90 * t125 + t114, t110; -sin(qJ(1)) * t122 + t99, t99, t100 * t90 - t126 * t131 - t109, t90 * t119 + t91 * t125, t109; 0, 0, -qJD(3) * t131 - t108 + t118, t120, t108;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end