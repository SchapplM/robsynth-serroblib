% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% JgD_rot [3x7]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S7RRRRRRR1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobigD_rot_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobigD_rot_sym_varpar: qJD has to be [7x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S7RRRRRRR1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobigD_rot_sym_varpar: pkin has to be [4x1] (double)');
%% Function calls
if link_index == 0
	JgD_rot=S7RRRRRRR1_jacobigD_rot_0_sym_varpar(qJ, qJD, pkin);
elseif link_index == 1
	JgD_rot=S7RRRRRRR1_jacobigD_rot_1_sym_varpar(qJ, qJD, pkin);
elseif link_index == 2
	JgD_rot=S7RRRRRRR1_jacobigD_rot_2_sym_varpar(qJ, qJD, pkin);
elseif link_index == 3
	JgD_rot=S7RRRRRRR1_jacobigD_rot_3_sym_varpar(qJ, qJD, pkin);
elseif link_index == 4
	JgD_rot=S7RRRRRRR1_jacobigD_rot_4_sym_varpar(qJ, qJD, pkin);
elseif link_index == 5
	JgD_rot=S7RRRRRRR1_jacobigD_rot_5_sym_varpar(qJ, qJD, pkin);
elseif link_index == 6
	JgD_rot=S7RRRRRRR1_jacobigD_rot_6_sym_varpar(qJ, qJD, pkin);
elseif link_index == 7
	JgD_rot=S7RRRRRRR1_jacobigD_rot_7_sym_varpar(qJ, qJD, pkin);
else
	JgD_rot=NaN(3,7);
end