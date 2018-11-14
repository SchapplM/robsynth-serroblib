% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% 
% Output:
% JR_rot [9x3]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JR_rot = S3RPP1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),uint8(0),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_jacobiR_rot_sym_varpar: qJ has to be [3x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S3RPP1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_jacobiR_rot_sym_varpar: pkin has to be [3x1] (double)');
%% Function calls
if link_index == 0
	JR_rot=S3RPP1_jacobiR_rot_0_sym_varpar(qJ, pkin);
elseif link_index == 1
	JR_rot=S3RPP1_jacobiR_rot_1_sym_varpar(qJ, pkin);
elseif link_index == 2
	JR_rot=S3RPP1_jacobiR_rot_2_sym_varpar(qJ, pkin);
elseif link_index == 3
	JR_rot=S3RPP1_jacobiR_rot_3_sym_varpar(qJ, pkin);
else
	JR_rot=NaN(9,3);
end