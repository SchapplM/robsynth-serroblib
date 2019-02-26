% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:48
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPPRRR2_jacobigD_rot_floatb_twist_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_jacobigD_rot_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_jacobigD_rot_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR2_jacobigD_rot_floatb_twist_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_jacobigD_rot_floatb_twist_sym_varpar: pkin has to be [11x1] (double)');
%% Function calls
if link_index == 0
	JgD_rot=S6RPPRRR2_jacobigD_rot_0_floatb_twist_sym_varpar(qJ, qJD, pkin);
elseif link_index == 1
	JgD_rot=S6RPPRRR2_jacobigD_rot_1_floatb_twist_sym_varpar(qJ, qJD, pkin);
elseif link_index == 2
	JgD_rot=S6RPPRRR2_jacobigD_rot_2_floatb_twist_sym_varpar(qJ, qJD, pkin);
elseif link_index == 3
	JgD_rot=S6RPPRRR2_jacobigD_rot_3_floatb_twist_sym_varpar(qJ, qJD, pkin);
elseif link_index == 4
	JgD_rot=S6RPPRRR2_jacobigD_rot_4_floatb_twist_sym_varpar(qJ, qJD, pkin);
elseif link_index == 5
	JgD_rot=S6RPPRRR2_jacobigD_rot_5_floatb_twist_sym_varpar(qJ, qJD, pkin);
elseif link_index == 6
	JgD_rot=S6RPPRRR2_jacobigD_rot_6_floatb_twist_sym_varpar(qJ, qJD, pkin);
else
	JgD_rot=NaN(6,6);
end