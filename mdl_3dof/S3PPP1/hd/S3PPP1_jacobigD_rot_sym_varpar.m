% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S3PPP1
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,theta1]';
% 
% Output:
% JgD_rot [3x3]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-17 09:48
% Revision: 3acd05283b8979b361f80d69cfa1c98d98241298 (2019-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S3PPP1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),uint8(0),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPP1_jacobigD_rot_sym_varpar: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPP1_jacobigD_rot_sym_varpar: qJD has to be [3x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S3PPP1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPP1_jacobigD_rot_sym_varpar: pkin has to be [3x1] (double)');
%% Function calls
if link_index == 0
	JgD_rot=S3PPP1_jacobigD_rot_0_sym_varpar(qJ, qJD, pkin);
elseif link_index == 1
	JgD_rot=S3PPP1_jacobigD_rot_1_sym_varpar(qJ, qJD, pkin);
elseif link_index == 2
	JgD_rot=S3PPP1_jacobigD_rot_2_sym_varpar(qJ, qJD, pkin);
elseif link_index == 3
	JgD_rot=S3PPP1_jacobigD_rot_3_sym_varpar(qJ, qJD, pkin);
else
	JgD_rot=NaN(3,3);
end