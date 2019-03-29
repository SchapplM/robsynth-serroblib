% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% JgD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S5RRRRR2_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobigD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_jacobigD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR2_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobigD_rot_sym_varpar: pkin has to be [2x1] (double)');
%% Function calls
if link_index == 0
	JgD_rot=S5RRRRR2_jacobigD_rot_0_sym_varpar(qJ, qJD, pkin);
elseif link_index == 1
	JgD_rot=S5RRRRR2_jacobigD_rot_1_sym_varpar(qJ, qJD, pkin);
elseif link_index == 2
	JgD_rot=S5RRRRR2_jacobigD_rot_2_sym_varpar(qJ, qJD, pkin);
elseif link_index == 3
	JgD_rot=S5RRRRR2_jacobigD_rot_3_sym_varpar(qJ, qJD, pkin);
elseif link_index == 4
	JgD_rot=S5RRRRR2_jacobigD_rot_4_sym_varpar(qJ, qJD, pkin);
elseif link_index == 5
	JgD_rot=S5RRRRR2_jacobigD_rot_5_sym_varpar(qJ, qJD, pkin);
else
	JgD_rot=NaN(3,5);
end