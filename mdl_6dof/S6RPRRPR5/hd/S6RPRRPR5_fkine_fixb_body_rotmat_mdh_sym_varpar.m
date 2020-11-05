% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:48
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:47:59
	% EndTime: 2020-11-04 21:47:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:47:59
	% EndTime: 2020-11-04 21:47:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t59 = cos(qJ(1));
	t58 = sin(qJ(1));
	t1 = [t59, -t58, 0, 0; t58, t59, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:47:59
	% EndTime: 2020-11-04 21:47:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t63 = cos(qJ(1));
	t62 = sin(qJ(1));
	t61 = cos(pkin(10));
	t60 = sin(pkin(10));
	t1 = [t63 * t61, -t63 * t60, t62, t63 * pkin(1) + t62 * qJ(2) + 0; t62 * t61, -t62 * t60, -t63, t62 * pkin(1) - t63 * qJ(2) + 0; t60, t61, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:47:59
	% EndTime: 2020-11-04 21:47:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t70 = cos(qJ(1));
	t69 = sin(qJ(1));
	t68 = pkin(7) + qJ(2);
	t67 = pkin(10) + qJ(3);
	t66 = cos(t67);
	t65 = sin(t67);
	t64 = cos(pkin(10)) * pkin(2) + pkin(1);
	t1 = [t70 * t66, -t70 * t65, t69, t70 * t64 + t68 * t69 + 0; t69 * t66, -t69 * t65, -t70, t69 * t64 - t70 * t68 + 0; t65, t66, 0, sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:47:59
	% EndTime: 2020-11-04 21:47:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->15), mult. (14->12), div. (0->0), fcn. (22->8), ass. (0->9)
	t76 = pkin(10) + qJ(3);
	t78 = cos(qJ(1));
	t77 = sin(qJ(1));
	t75 = -pkin(8) - pkin(7) - qJ(2);
	t74 = qJ(4) + t76;
	t73 = cos(t74);
	t72 = sin(t74);
	t71 = pkin(3) * cos(t76) + cos(pkin(10)) * pkin(2) + pkin(1);
	t1 = [t78 * t73, -t78 * t72, t77, t78 * t71 - t77 * t75 + 0; t77 * t73, -t77 * t72, -t78, t77 * t71 + t78 * t75 + 0; t72, t73, 0, pkin(3) * sin(t76) + sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:47:59
	% EndTime: 2020-11-04 21:47:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (53->21), mult. (24->16), div. (0->0), fcn. (32->8), ass. (0->9)
	t84 = pkin(10) + qJ(3);
	t82 = qJ(4) + t84;
	t80 = sin(t82);
	t81 = cos(t82);
	t87 = pkin(4) * t81 + qJ(5) * t80 + pkin(3) * cos(t84) + cos(pkin(10)) * pkin(2) + pkin(1);
	t86 = cos(qJ(1));
	t85 = sin(qJ(1));
	t83 = -pkin(8) - pkin(7) - qJ(2);
	t1 = [t85, -t86 * t81, t86 * t80, -t85 * t83 + t87 * t86 + 0; -t86, -t85 * t81, t85 * t80, t86 * t83 + t87 * t85 + 0; 0, -t80, -t81, t80 * pkin(4) - t81 * qJ(5) + pkin(3) * sin(t84) + sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:47:59
	% EndTime: 2020-11-04 21:47:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (65->23), mult. (38->24), div. (0->0), fcn. (51->10), ass. (0->16)
	t95 = sin(qJ(6));
	t96 = sin(qJ(1));
	t104 = t96 * t95;
	t97 = cos(qJ(6));
	t103 = t96 * t97;
	t98 = cos(qJ(1));
	t102 = t98 * t95;
	t101 = t98 * t97;
	t94 = pkin(10) + qJ(3);
	t92 = qJ(4) + t94;
	t88 = sin(t92);
	t89 = cos(t92);
	t99 = pkin(4) + pkin(9);
	t100 = pkin(3) * cos(t94) + qJ(5) * t88 + t89 * t99 + cos(pkin(10)) * pkin(2) + pkin(1);
	t93 = qJ(2) + pkin(5) + pkin(7) + pkin(8);
	t1 = [t88 * t102 + t103, t88 * t101 - t104, t98 * t89, t100 * t98 + t93 * t96 + 0; t88 * t104 - t101, t88 * t103 + t102, t96 * t89, t100 * t96 - t93 * t98 + 0; -t89 * t95, -t89 * t97, t88, t99 * t88 - t89 * qJ(5) + sin(pkin(10)) * pkin(2) + pkin(3) * sin(t94) + 0 + pkin(6); 0, 0, 0, 1;];
	Tc_mdh = t1;
end