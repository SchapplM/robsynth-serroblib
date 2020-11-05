% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPPR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:35
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPPR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:35:06
	% EndTime: 2020-11-04 21:35:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:35:06
	% EndTime: 2020-11-04 21:35:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t1 = [t56, -t55, 0, 0; t55, t56, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:35:06
	% EndTime: 2020-11-04 21:35:06
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t58 = cos(qJ(1));
	t57 = sin(qJ(1));
	t1 = [0, -t58, t57, t58 * pkin(1) + t57 * qJ(2) + 0; 0, -t57, -t58, t57 * pkin(1) - t58 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:35:06
	% EndTime: 2020-11-04 21:35:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t63 = pkin(1) + pkin(7);
	t62 = cos(qJ(1));
	t61 = cos(qJ(3));
	t60 = sin(qJ(1));
	t59 = sin(qJ(3));
	t1 = [t60 * t59, t60 * t61, t62, t60 * qJ(2) + t63 * t62 + 0; -t62 * t59, -t62 * t61, t60, -t62 * qJ(2) + t63 * t60 + 0; t61, -t59, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:35:06
	% EndTime: 2020-11-04 21:35:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->13), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t70 = cos(qJ(1));
	t69 = sin(qJ(1));
	t68 = qJ(3) + pkin(9);
	t67 = pkin(1) + pkin(7) + qJ(4);
	t66 = cos(t68);
	t65 = sin(t68);
	t64 = sin(qJ(3)) * pkin(3) + qJ(2);
	t1 = [t69 * t65, t69 * t66, t70, t64 * t69 + t67 * t70 + 0; -t70 * t65, -t70 * t66, t69, -t64 * t70 + t67 * t69 + 0; t66, -t65, 0, cos(qJ(3)) * pkin(3) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:35:06
	% EndTime: 2020-11-04 21:35:06
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (34->19), mult. (23->17), div. (0->0), fcn. (31->8), ass. (0->11)
	t77 = sin(pkin(9));
	t78 = cos(pkin(9));
	t81 = cos(qJ(3));
	t83 = (pkin(4) * t78 + qJ(5) * t77 + pkin(3)) * sin(qJ(3)) - (-t77 * pkin(4) + qJ(5) * t78) * t81 + qJ(2);
	t82 = cos(qJ(1));
	t80 = sin(qJ(1));
	t76 = qJ(3) + pkin(9);
	t75 = pkin(1) + pkin(7) + qJ(4);
	t74 = cos(t76);
	t73 = sin(t76);
	t1 = [t82, -t80 * t73, -t80 * t74, t75 * t82 + t83 * t80 + 0; t80, t82 * t73, t82 * t74, t75 * t80 - t83 * t82 + 0; 0, -t74, t73, t81 * pkin(3) + t74 * pkin(4) + t73 * qJ(5) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:35:06
	% EndTime: 2020-11-04 21:35:06
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (46->22), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->18)
	t89 = sin(pkin(9));
	t90 = cos(pkin(9));
	t95 = cos(qJ(3));
	t97 = pkin(4) + pkin(8);
	t104 = (qJ(5) * t89 + t97 * t90 + pkin(3)) * sin(qJ(3)) - (qJ(5) * t90 - t89 * t97) * t95 + qJ(2);
	t91 = sin(qJ(6));
	t93 = sin(qJ(1));
	t103 = t93 * t91;
	t94 = cos(qJ(6));
	t102 = t93 * t94;
	t96 = cos(qJ(1));
	t101 = t96 * t91;
	t100 = t96 * t94;
	t88 = qJ(3) + pkin(9);
	t87 = qJ(4) + pkin(1) + pkin(5) + pkin(7);
	t86 = cos(t88);
	t85 = sin(t88);
	t1 = [-t86 * t103 + t100, -t86 * t102 - t101, t93 * t85, t104 * t93 + t87 * t96 + 0; t86 * t101 + t102, t86 * t100 - t103, -t96 * t85, -t104 * t96 + t87 * t93 + 0; t85 * t91, t85 * t94, t86, t95 * pkin(3) + t85 * qJ(5) + t97 * t86 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end