% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:48
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:48:40
	% EndTime: 2020-11-04 21:48:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:48:40
	% EndTime: 2020-11-04 21:48:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t55 = cos(qJ(1));
	t54 = sin(qJ(1));
	t1 = [t55, -t54, 0, 0; t54, t55, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:48:40
	% EndTime: 2020-11-04 21:48:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t1 = [0, -t57, t56, t57 * pkin(1) + t56 * qJ(2) + 0; 0, -t56, -t57, t56 * pkin(1) - t57 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:48:40
	% EndTime: 2020-11-04 21:48:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t62 = pkin(1) + pkin(7);
	t61 = cos(qJ(1));
	t60 = cos(qJ(3));
	t59 = sin(qJ(1));
	t58 = sin(qJ(3));
	t1 = [t59 * t58, t59 * t60, t61, t59 * qJ(2) + t62 * t61 + 0; -t61 * t58, -t61 * t60, t59, -t61 * qJ(2) + t62 * t59 + 0; t60, -t58, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:48:40
	% EndTime: 2020-11-04 21:48:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->13), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t69 = cos(qJ(1));
	t68 = sin(qJ(1));
	t67 = qJ(3) + qJ(4);
	t66 = pkin(1) + pkin(7) + pkin(8);
	t65 = cos(t67);
	t64 = sin(t67);
	t63 = sin(qJ(3)) * pkin(3) + qJ(2);
	t1 = [t68 * t64, t68 * t65, t69, t63 * t68 + t66 * t69 + 0; -t69 * t64, -t69 * t65, t68, -t63 * t69 + t66 * t68 + 0; t65, -t64, 0, cos(qJ(3)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:48:41
	% EndTime: 2020-11-04 21:48:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (36->16), mult. (16->12), div. (0->0), fcn. (24->8), ass. (0->9)
	t76 = qJ(3) + qJ(4);
	t79 = pkin(4) * sin(t76) + sin(qJ(3)) * pkin(3) + qJ(2);
	t78 = cos(qJ(1));
	t77 = sin(qJ(1));
	t75 = qJ(5) + pkin(1) + pkin(7) + pkin(8);
	t73 = pkin(10) + t76;
	t71 = cos(t73);
	t70 = sin(t73);
	t1 = [t77 * t70, t77 * t71, t78, t75 * t78 + t79 * t77 + 0; -t78 * t70, -t78 * t71, t77, t75 * t77 - t79 * t78 + 0; t71, -t70, 0, pkin(4) * cos(t76) + cos(qJ(3)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:48:41
	% EndTime: 2020-11-04 21:48:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (63->23), mult. (38->24), div. (0->0), fcn. (51->10), ass. (0->15)
	t87 = sin(qJ(6));
	t88 = sin(qJ(1));
	t95 = t88 * t87;
	t89 = cos(qJ(6));
	t94 = t88 * t89;
	t90 = cos(qJ(1));
	t93 = t90 * t87;
	t92 = t90 * t89;
	t86 = qJ(3) + qJ(4);
	t83 = pkin(10) + t86;
	t80 = sin(t83);
	t81 = cos(t83);
	t91 = pkin(4) * sin(t86) + pkin(5) * t80 - pkin(9) * t81 + sin(qJ(3)) * pkin(3) + qJ(2);
	t85 = qJ(5) + pkin(1) + pkin(7) + pkin(8);
	t1 = [t80 * t94 + t93, -t80 * t95 + t92, -t88 * t81, t85 * t90 + t91 * t88 + 0; -t80 * t92 + t95, t80 * t93 + t94, t90 * t81, t85 * t88 - t91 * t90 + 0; t81 * t89, -t81 * t87, t80, t81 * pkin(5) + t80 * pkin(9) + pkin(4) * cos(t86) + cos(qJ(3)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end