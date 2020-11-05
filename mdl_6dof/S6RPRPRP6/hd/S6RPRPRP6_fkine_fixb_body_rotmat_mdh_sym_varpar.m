% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRP6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:37
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:29
	% EndTime: 2020-11-04 21:37:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:29
	% EndTime: 2020-11-04 21:37:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t61 = cos(qJ(1));
	t60 = sin(qJ(1));
	t1 = [t61, -t60, 0, 0; t60, t61, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:29
	% EndTime: 2020-11-04 21:37:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t65 = cos(qJ(1));
	t64 = sin(qJ(1));
	t63 = cos(pkin(9));
	t62 = sin(pkin(9));
	t1 = [t65 * t63, -t65 * t62, t64, t65 * pkin(1) + t64 * qJ(2) + 0; t64 * t63, -t64 * t62, -t65, t64 * pkin(1) - t65 * qJ(2) + 0; t62, t63, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:29
	% EndTime: 2020-11-04 21:37:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t72 = cos(qJ(1));
	t71 = sin(qJ(1));
	t70 = pkin(7) + qJ(2);
	t69 = pkin(9) + qJ(3);
	t68 = cos(t69);
	t67 = sin(t69);
	t66 = cos(pkin(9)) * pkin(2) + pkin(1);
	t1 = [t72 * t68, -t72 * t67, t71, t72 * t66 + t70 * t71 + 0; t71 * t68, -t71 * t67, -t72, t71 * t66 - t72 * t70 + 0; t67, t68, 0, sin(pkin(9)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:29
	% EndTime: 2020-11-04 21:37:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->18), mult. (21->14), div. (0->0), fcn. (29->6), ass. (0->8)
	t76 = pkin(9) + qJ(3);
	t74 = sin(t76);
	t75 = cos(t76);
	t80 = pkin(3) * t75 + qJ(4) * t74 + cos(pkin(9)) * pkin(2) + pkin(1);
	t79 = cos(qJ(1));
	t78 = sin(qJ(1));
	t77 = pkin(7) + qJ(2);
	t1 = [t78, -t79 * t75, t79 * t74, t77 * t78 + t80 * t79 + 0; -t79, -t78 * t75, t78 * t74, -t79 * t77 + t80 * t78 + 0; 0, -t74, -t75, t74 * pkin(3) - t75 * qJ(4) + sin(pkin(9)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:29
	% EndTime: 2020-11-04 21:37:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (44->22), mult. (37->26), div. (0->0), fcn. (50->10), ass. (0->17)
	t88 = sin(qJ(5));
	t89 = sin(qJ(1));
	t96 = t89 * t88;
	t90 = cos(qJ(5));
	t95 = t89 * t90;
	t91 = cos(qJ(1));
	t94 = t91 * t88;
	t93 = t91 * t90;
	t92 = pkin(3) + pkin(8);
	t87 = cos(pkin(9));
	t86 = sin(pkin(9));
	t85 = pkin(9) + qJ(3);
	t84 = qJ(2) + pkin(4) + pkin(7);
	t83 = cos(t85);
	t82 = sin(t85);
	t81 = (qJ(4) * t86 + t92 * t87) * cos(qJ(3)) + (qJ(4) * t87 - t86 * t92) * sin(qJ(3)) + t87 * pkin(2) + pkin(1);
	t1 = [t82 * t94 + t95, t82 * t93 - t96, t91 * t83, t81 * t91 + t84 * t89 + 0; t82 * t96 - t93, t82 * t95 + t94, t89 * t83, t81 * t89 - t84 * t91 + 0; -t83 * t88, -t83 * t90, t82, t86 * pkin(2) - t83 * qJ(4) + t92 * t82 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:29
	% EndTime: 2020-11-04 21:37:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (56->23), mult. (50->24), div. (0->0), fcn. (63->8), ass. (0->16)
	t115 = pkin(3) + qJ(6) + pkin(8);
	t106 = cos(qJ(5));
	t114 = t106 * pkin(5) + pkin(4) + pkin(7) + qJ(2);
	t104 = sin(qJ(5));
	t105 = sin(qJ(1));
	t113 = t105 * t104;
	t112 = t105 * t106;
	t107 = cos(qJ(1));
	t111 = t107 * t104;
	t110 = t107 * t106;
	t109 = pkin(5) * t104 + qJ(4);
	t101 = pkin(9) + qJ(3);
	t100 = cos(t101);
	t99 = sin(t101);
	t108 = t115 * t100 + t109 * t99 + cos(pkin(9)) * pkin(2) + pkin(1);
	t1 = [t99 * t111 + t112, t99 * t110 - t113, t107 * t100, t114 * t105 + t108 * t107 + 0; t99 * t113 - t110, t99 * t112 + t111, t105 * t100, t108 * t105 - t114 * t107 + 0; -t100 * t104, -t100 * t106, t99, sin(pkin(9)) * pkin(2) + pkin(6) + 0 + t115 * t99 - t109 * t100; 0, 0, 0, 1;];
	Tc_mdh = t1;
end